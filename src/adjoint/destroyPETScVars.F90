!
!     ******************************************************************
!     *                                                                *
!     * File:          destroyPETScVars.F90                            *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 08-15-2005                                      *
!     * Last modified: 04-08-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine destroyPETScVars
!
!     ******************************************************************
!     *                                                                *
!     * Free work space. All PETSc objects should be destroyed when    *
!     * they are no longer needed.                                     *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use communication
      use flowVarRefState ! magnetic
      implicit none
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! Write a message that the PETSc objects are being destroyed.
      ! Of course only processor 0 does this.

      if(myID == 0) then
        print "(a)", "#"
        print "(a)", "# Destroying PETSc objects ..."
      endif

      ! VecDestroy - Destroys a vector.
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! call VecDestroy(Vec v, PetscErrorCode ierr)
      !
      ! Collective on Vec
      !
      ! Input Parameters
      !   v -the vector
      !
      ! see .../petsc/docs/manualpages/Vec/VecDestroy.html
      ! or PETSc users manual, pp.37

      call VecDestroy(psi, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector psi")

      call VecDestroy(dJdW, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector dJdW")

      call VecDestroy(pvr, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector pvr")

      call VecDestroy(dJda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector dJda")

      call VecDestroy(dIda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector dIda")

      call VecDestroy(dJdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector dJdx")

      call VecDestroy(dIdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector dIdx")

      ! MatDestroy - Frees space taken by a matrix.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatDestroy(Mat A, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameter
      !   A -the matrix 
      !
      ! see .../petsc/docs/manualpages/Mat/MatDestroy.html
      ! or PETSc users manual, pp.61

      call MatDestroy(dRdW, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRdW")

      call MatDestroy(dRda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRda")

      call MatDestroy(dRdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRdx")

      ! KSPDestroy - Destroys KSP context.
      !
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPDestroy(KSP ksp, PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameter
      !   ksp -iterative context obtained from KSPCreate() 
      !
      ! see .../petsc/docs/manualpages/KSP/KSPDestroy.html
      ! or PETSc users manual, pp.64

      call KSPDestroy(ksp, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy KSP context")

      ! Synchronize the processors.

      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

#endif

      end subroutine destroyPETScVars
