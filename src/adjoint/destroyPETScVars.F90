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
      use flowVarRefState ! 
      use inputADjoint    !ApproxPC
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

      !print *,'vecs'

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

      call VecDestroy(dJcdW, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector dJcdW")

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

      call VecDestroy(dIdxsDV, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector dIdxsDV")

      call VecDestroy(phic, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector phic")

      call VecDestroy(dJdc, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector djdc")
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
      !print *,'mats'

!!$      call MatDestroy(dRdW, PETScIerr)
!!$
!!$      if( PETScIerr/=0 ) &
!!$        call terminate("destroyPETScVars", &
!!$                       "Could not destroy matrix dRdW")

      call MatDestroy(dRdWt, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRdWt")
!!$
!!$      call MatDestroy(dRdWPre, PETScIerr)
!!$
!!$      if( PETScIerr/=0 ) &
!!$        call terminate("destroyPETScVars", &
!!$                       "Could not destroy matrix dRdWpre")
      if (ApproxPC) then
         !print *,'mats,prpwpre'
         call MatDestroy(dRdWpret, PETScIerr)
         
         if( PETScIerr/=0 ) &
              call terminate("destroyPETScVars", &
              "Could not destroy matrix dRdWpret")
      endif
      !print *,'mats,drda'
      call MatDestroy(dRda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRda")
      !print *,'mats,drdx'
      call MatDestroy(dRdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRdx")
      !print *,'mats,dcdw'
      call MatDestroy(dCdw, PETScIerr)

      if( PETScIerr/=0 ) &
           call terminate("destroyPETScVars", &
           "Could not destroy matrix dcdw")
      !print *,'matsdcdx'
      call MatDestroy(dCdx, PETScIerr)

      if( PETScIerr/=0 ) &
           call terminate("destroyPETScVars", &
           "Could not destroy matrix dcdx")
      !print *,'matsdcda'
      call MatDestroy(dCda, PETScIerr)

      if( PETScIerr/=0 ) &
           call terminate("destroyPETScVars", &
           "Could not destroy matrix dcda")
      !print *,'mats finished'

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
      !print *,'ksp'
      call KSPDestroy(ksp, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy KSP context")

      ! Synchronize the processors.

      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

#endif

      end subroutine destroyPETScVars
