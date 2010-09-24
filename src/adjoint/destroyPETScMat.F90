!
!     ******************************************************************
!     *                                                                *
!     * File:          destroyPETScMat.F90                             *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 05-13-2010                                      *
!     * Last modified: 05-13-2010                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine destroyPETScMat
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
        print "(a)", "# Destroying PETSc MAT objects ..."
      endif

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

      call MatDestroy(dRdWt, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRdWt")

      call MatDestroy(dRdWPre, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRdWpre")

      call MatDestroy(dRdWpret, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRdWpret")


      call MatDestroy(dRda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRda")

      call MatDestroy(dRdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRdx")

      call MatDestroy(dCdw, PETScIerr)

      if( PETScIerr/=0 ) &
           call terminate("destroyPETScVars", &
           "Could not destroy matrix dcdw")

      call MatDestroy(dCdx, PETScIerr)

      if( PETScIerr/=0 ) &
           call terminate("destroyPETScVars", &
           "Could not destroy matrix dcdx")

      call MatDestroy(dCda, PETScIerr)

      if( PETScIerr/=0 ) &
           call terminate("destroyPETScVars", &
           "Could not destroy matrix dcda")


      ! Synchronize the processors.

      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

#endif

    end subroutine destroyPETScMat
