!
!     ******************************************************************
!     *                                                                *
!     * File:          computeADjointGradientSpatial.F90               *
!     * Author:        Andre C. Marta, C.A.(Sandy) Mader               *
!     * Starting date: 11-29-2006                                      *
!     * Last modified: 05-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine computeADjointGradientSpatial(costFunction)
!
!     ******************************************************************
!     *                                                                *
!     * Computes the total gradient of the cost function J with        *
!     * respect to the spatial design variables using the adjoint sol. *
!     * In addition, it also copies the gradient from the PETSc object *
!     * to the storage array functionGrad in the root processor.       *
!     *                                                                *
!     * Compute the total cost/constraint function sensitivity as      *
!     * dIdx = dJdx - psi^T dRdx, where                                *
!     *                                                                *
!     *  dIdx [3*nNodesGlobal] is the total gradient of the cost       *
!     *            function with respect to the design variables,      *
!     *                                                                *
!     *  dJdx [3*nNodesglobal] is the partial gradient of the cost     *
!     *            function with respect to the design variables,      *
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
#ifndef USE_NO_PETSC

      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Computing total sensitivity wrt Spatial..."

      ! Get the initial time.

      call cpu_time(time(1))

      ! Compute the second contribution term psi^T dRdx, which requires
      ! a matrix-vector multiplication.

      call MatMultTranspose(dRdx,psi,dIdx,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientSpatial", &
                       "Error in MatMultTranspose X")

      ! Add the first dJ/dx term to get the total gradient and take into
      ! account the negative sign in front of the second contribution.

      call VecAYPX(dIdx,PETScNegOne,dJdx,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientSpatial", &
                       "Error in VecAYPX X")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Computing total sensitivity wrt Spatial time (s) =", &
                    timeAdj

      ! View the solution vector dIdx.
 
     if(  debug ) then

        if( PETScRank==0 ) then
          write(*,*) "# ============================ "
          write(*,*) "#  dIdx total gradient vector  "
          write(*,*) "# ============================ "
        endif

        call VecView(dIdx,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
!        call VecView(dIdx,PETSC_VIEWER_DRAW_WORLD,PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("computeADjointGradientSpatial", &
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
	
   !   print *,'irange',ilow,ihigh,petscrank
   !   print *,'irange',iLow,iHigh

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientSpatial", &
                       "Error in VecGetOwnershipRange dJdx")
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

        call VecGetValues(dIdx, 1, idxmg, &
                          functionGradLocal(n), PETScIerr)
	
	!functionGradLocal(n) = idxmg !to debug interface with meshwarping derivatives	

        if( PETScIerr/=0 ) then
          write(errorMessage,99) &
                "Error in VecGetValues for global node", idxmg
          call terminate("computeADjointGradientSpatial", errorMessage)
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

      !print *,'costfunction',costFunction,nDisplsGlobal,'shapes',shape(functionGradSpatial),shape(functionGradLocal)

!      call mpi_gatherv(functionGradLocal, nDesignLocal, sumb_real, &
!                       functionGradSpatial(costFunction,:), nDesignGlobal,&
!                       nDisplsGlobal, sumb_real, &
!                       0, PETSC_COMM_WORLD, PETScIerr)
       call mpi_allgatherv(functionGradLocal, nDesignLocal, sumb_real, &
                       functionGradSpatial(costFunction,:), nDesignGlobal,&
                       nDisplsGlobal, sumb_real, &
                        PETSC_COMM_WORLD, PETScIerr)
!      globalNode => flowdoms(3,1,1)%globalNode
!      idx = globalnode(2,7,13)*3+3 

!<<<<<<< .mine
!
!      globalNode =>flowdoms(7,1,1)%globalNode
!!      globalNode =>flowdoms(3,1,1)%globalNode
!!###      globalNode =>flowdoms(1,1,1)%globalNode
!=======
!********************************************************
      !use to print out values to chack accuracy in multiblock cases 
!      globalNode =>flowdoms(1,1,1)%globalNode 
      !idxlocal = globalNode(4,10,10)*3+3 !Obliquewing coarse block 7
!      idxlocal = globalNode(3,9,3)*3+3 !Obliquewing 20k coarse block 30
      !print *, 'globalnode',idxlocal,petscrank
      
      call mpi_gather(idxLocal, 1, sumb_integer, &
                      idxGlobal, 1, sumb_integer, &
                      0, PETSC_COMM_WORLD, PETScIerr)

     ! if(PetscRank == 0)then
!	 do i =1,petscsize
!         print *,'costfunction derivative',idxglobal(i),functionGradSpatial(costFunction,idxglobal(i)),functionGradSpatial(costFunction,46485),i
!	 enddo
!      endif
!**********************************************************

!***********  
!original single block verification code.   
!***************
!      globalNode =>flowdoms(7,1,1)%globalNode !Obliquewing coarse block 7
!      globalNode =>flowdoms(8,1,1)%globalNode !Obliquewing coarse block 8
!      globalNode =>flowdoms(1,1,1)%globalNode !infinite wing
!>>>>>>> .r625
!      idx = globalNode(2,7,13)*3+3 !(location times 3 for xyzdof + 3 for current node
!      idx = globalNode(10,1,3)*3+3 !infinite wing block 1
!      idx = globalNode(4,10,10)*3+3 !Obliquewing coarse block 7
!      idx = globalNode(4,1,10)*3+3 !Obliquewing coarse block 8

!      print *,'costfunction derivative',idx,functionGradSpatial(costFunction,idx)
!****************


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

      end subroutine computeADjointGradientSpatial
