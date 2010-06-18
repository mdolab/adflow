
!
!     ******************************************************************
!     *                                                                *
!     * File:          setupADjointRHSStability.F90                    *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 11-27-2009                                      *
!     * Last modified: 11-27-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupADjointRHSStability(level,costFunction)

!
!     ******************************************************************
!     *                                                                *
!     * Compute the right hand side of the discrete ADjoint problem    *
!     * in question. Notice that this right hand side is problem /     *
!     * cost function J dependent.                                     *
!     *                                                                *
!     * The ordering of the unknowns in the ADjoint vector used here   *
!     * is based on the global node numbering and is consistent with   *
!     * the ordering used in the matrix for the ADjoint problem        *
!     * assembled in setupADjointMatrix.                               *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level, costFunction
!
!     Local variables.
!

      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj
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
        write(*,10) "Assembling ADjoint Stability RHS vector..."

      ! Get the initial time.

      call cpu_time(time(1))

      ! Reset the RHS vector dJ/dW by assigning the value zero to all
      ! its components.

      ! VecSet - Sets all components of vector to a single scalar value.
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! call VecSet(Vec x,PetscScalar alpha, PetscErrorCode ierr)
      !
      ! Collective on Vec
      !
      ! Input Parameters
      !   x     - the vector
      !   alpha	- the scalar
      !
      ! Output Parameter
      !   x -the vector
      !
      ! Note
      ! For a vector of dimension n, VecSet() computes
      ! x[i] = alpha, for i=1,...,n,
      ! so that all vector entries then equal the identical scalar
      ! value, alpha. Use the more general routine VecSetValues() to
      ! set different vector entries.
      !
      ! You CANNOT call this after you have called VecSetValues() but
      ! before you call VecAssemblyBegin/End(). 
      !
      ! see .../petsc/docs/manualpages/Vec/VecSet.html
      ! or PETSc users manual, pp.36

      call VecSet(dJdC,PETScZero,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointRHSStability", "Error in VecSet")

      !zero the matrix for dCdW Insert call
      call MatZeroEntries(dCdw,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointRHSStability",&
        "Error in MatZeroEntries dCdw")

      !Compute dCdw 
      call setupADjointdCdwStability(level,costFunction)

      !compute dIdc
      call setupADjointdIdCStability(level,costFunction)

      !multiply to get djdw
      call MatMultTranspose(dCdw,dJdc,dJdw,PETScIerr)
    

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, SUMB_PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling ADjoint TS RHS vector time (s) = ", timeAdj


      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)

#endif

    end subroutine setupADjointRHSStability
