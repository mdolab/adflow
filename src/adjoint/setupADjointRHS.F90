
!
!     ******************************************************************
!     *                                                                *
!     * File:          setupADjointRHS.F90                             *
!     * Author:        Andre C. Marta, C.A.(Sandy) Mader               *
!     * Starting date: 07-27-2006                                      *
!     * Last modified: 05-13-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupADjointRHS(level,sps,costFunction)
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
      integer(kind=intType), intent(in) :: level, sps, costFunction
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
        write(*,10) "Assembling ADjoint RHS vector..."

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

      call VecSet(dJdW,PETScZero,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointRHS", "Error in VecSet")
!
!     ******************************************************************
!     *                                                                *
!     * Select case over different choices of cost functions.          *
!     * The following routine calls fill in the PETSc vector dJdW by   *
!     * making calls to the PETSc routine VecSetValuesBlocked.         *
!     *                                                                *
!     ******************************************************************
!
      select case (costFunction)

        case (costFuncLiftCoef, &
              costFuncDragCoef, &
              costFuncForceXCoef,&
              costFuncForceYCoef,&
              costFuncForceZCoef,&
              costFuncMomXCoef, &
              costFuncMomYCoef, &
              costFuncMomZCoef)

          call setupADjointRHSAeroCoeff(level,costFunction)
!	call setupADjointRHSAeroCoeff(level,sps,costFunction)

        case default
          write(*,*) "Invalid cost function ", costFunction
          stop

      end select
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

      call VecAssemblyBegin(dJdW,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupASjointRHS", "Error in VecAssemblyBegin")

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

      call VecAssemblyEnd  (dJdW,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointRHS", "Error in VecAssemblyEnd")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling ADjoint RHS vector time (s) = ", timeAdj
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
	!call VecView(dJdW,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
	call VecView(dJdW,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("setupADjointRHS", "Error in VecView")
        pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)

#endif

    end subroutine setupADjointRHS
