
!
!     ******************************************************************
!     *                                                                *
!     * File:          setupADjointTotalStruct.f90                       *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 08-19-2008                                      *
!     * Last modified: 08-19-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine setupADjointTotalStruct(lenadjoint,structAdjoint,costfunction)
!
!     ******************************************************************
!     *                                                                *
!     * Compute the right hand side Augmentation of the adjoint problem*
!     * contributed by the structures.                                 *
!     *                                                                *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars
      use mdData              !mdNSurfNodes
      use communication       ! procHalo(currentLevel)%nProcSend, myID

      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType),intent(in)::lenadjoint
      integer(kind=intType) :: modFamID
      real(kind=realType), intent(in),dimension(3,lenadjoint) :: structAdjoint
      integer(kind=intType), intent(in) :: costFunction
            
!      real(kind=realType), intent(in),dimension(3,nSurfNodesLocal) :: structAdjoint
!
!     Local variables.
!
      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj

      integer(kind=intType)::n,m,idxsurf
!      integer(kind=intType) :: idx!tmp for fd

      character(len=2*maxStringLen) :: errorMessage

      integer(kind=intType) :: idxmg, iLow, iHigh
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
!#ifndef USE_NO_PETSC
!      !determine the number of surface nodes for coupling matrix
!      call mdCreateNSurfNodesLocal
      modFamID = max(0, 1_intType)


      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Assembling ADjoint Struct Totals vector..."

      ! Get the initial time.

      call cpu_time(time(1))


      do m=1,nSurfNodesLocal
         do n=1,3
                 
            idxSurf = (m-1)*3+n + (mdNsurfNodes(myID,modFamID)*3)

            if (structAdjoint(n,m).ne.0.0)then
               call VecSetValues(phic, 1, idxSurf-1,  &
                    structAdjoint(n,m), INSERT_VALUES, PETScIerr)
               if( PETScIerr/=0 ) &
                    print *,'matrix setting error'
               
            endif
         enddo
      enddo
!
!     ******************************************************************
!     *                                                                *
!     * Complete the PETSc vector assembly process.                    *
!     *                                                                *
!     ******************************************************************
!
      ! VecAssemblyBegin - Begins assembling the vector. This routine
      ! should be called after completing all calls to VecSetValues().
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! call VecAssemblyBegin(Vec vec, PetscErrorCode ierr)
      !
      ! Collective on Vec
      !
      ! Input Parameter
      !   vec -the vector 
      !
      ! see .../petsc/docs/manualpages/Vec/VecAssemblyBegin.html

      call VecAssemblyBegin(phic,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointTotalStruct", "Error in VecAssemblyBegin")

      ! VecAssemblyEnd - Completes assembling the vector. This routine
      ! should be called after VecAssemblyBegin().
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! call VecAssemblyEnd(Vec vec, PetscErrorCode ierr)
      !
      ! Collective on Vec
      !
      ! Input Parameter
      !   vec -the vector 
      !
      ! see .../petsc/docs/manualpages/Vec/VecAssemblyEnd.html

      call VecAssemblyEnd  (phic,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointTotalStruct", "Error in VecAssemblyEnd")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling ADjoint Total Struct vector time (s) = ", timeAdj
!
!     ******************************************************************
!     *                                                                *
!     * Visualize the assembled vector.                                *
!     *                                                                *
!     ******************************************************************
!
      ! VecView - Views a vector object.
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! PetscErrorCode PETSCVEC_DLLEXPORT VecView(Vec vec, &
      !                                              PetscViewer viewer)
      !
      ! Collective on Vec
      !
      ! Input Parameters
      !   v      - the vector
      !   viewer - an optional visualization context
      !
      ! Notes
      ! The available visualization contexts include
      !   PETSC_VIEWER_STDOUT_SELF  - standard output (default)
      !   PETSC_VIEWER_STDOUT_WORLD - synchronized standard output where
      !    only the first processor opens the file. All other processors
      !    send their data to the first processor to print.
      !
      ! see .../petsc/docs/manualpages/Vec/VecView.html
      ! or PETSc users manual, pp.36,148

      if( debug ) then
	call VecView(phic,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
	!call VecView(phic,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("setupADjointTotalStruct", "Error in VecView")
        pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      !now multiply with the dSdw derivatives
      
      ! Compute the second contribution term phic^T dSdw, which requires
      ! a matrix-vector multiplication.

      call MatMultTranspose(dSdx,phic,dJcdx,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("SetupADjointTotalStruct", &
                       "Error in MatMultTranspose X")

!
!     ******************************************************************
!     *                                                                *
!     * Transfer solution from PETSc context.                          *
!     *                                                                *
!     ******************************************************************
!

      ! Query about the ownership range.
      ! iHigh is one more than the last element stored locally.

      call VecGetOwnershipRange(dJcdx, iLow, iHigh, PETScIerr)
	
      if( PETScIerr/=0 ) &
        call terminate("setupADjointTotalStruct", &
                       "Error in VecGetOwnershipRange dJcdx")

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

        call VecGetValues(dJcdx, 1, idxmg, &
                          functionGradLocal(n), PETScIerr)

        if( PETScIerr/=0 ) then
          write(errorMessage,99) &
                "Error in VecGetValues for global node", idxmg
          call terminate("setupADjointTotalStruct", errorMessage)
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
                       functionGradStruct(costFunction,:), nDesignGlobal,&
                       nDisplsGlobal, sumb_real, &
                        PETSC_COMM_WORLD, PETScIerr)

      ! Release memory to store the local function gradient values.

      if (allocated(functionGradLocal)) deallocate(functionGradLocal)

      if (allocated(nDesignGlobal)) deallocate(nDesignGlobal)
      if (allocated(nDisplsGlobal)) deallocate(nDisplsGlobal)

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)
      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)
!#endif

    end subroutine setupADjointTotalStruct
