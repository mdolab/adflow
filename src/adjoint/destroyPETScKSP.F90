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
      subroutine destroyPETScKSP
!
!     ******************************************************************
!     *                                                                *
!     * Free KSP work space. All PETSc objects should be destroyed when*
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
        print "(a)", "# Destroying PETSc KSP Context ..."
      endif

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

      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

#endif

      end subroutine destroyPETScKSP
