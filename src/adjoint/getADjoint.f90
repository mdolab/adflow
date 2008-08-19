!
!     ******************************************************************
!     *                                                                *
!     * File:          getADjoint.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 08-18-2008                                      *
!     * Last modified: 08-18-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine getADjoint(costFunction)
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

      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Retrieving ADjoint Vector..."	
!
!     ******************************************************************
!     *                                                                *
!     * Transfer solution from PETSc context.                          *
!     *                                                                *
!     ******************************************************************
!

      ! Query about the ownership range.
      ! iHigh is one more than the last element stored locally.

      call VecGetOwnershipRange(psi, iLow, iHigh, PETScIerr)
	
!      print *,'irange',ilow,ihigh,petscrank
!      print *,'irange',iLow,iHigh

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

        call VecGetValues(psi, 1, idxmg, &
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


       call mpi_allgatherv(functionGradLocal, nDesignLocal, sumb_real, &
                       ADjoint(costFunction,:), nDesignGlobal,&
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

!#endif

    end subroutine getADjoint
