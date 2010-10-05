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

      !call VecDestroy(phic, PETScIerr)

!       if( PETScIerr/=0 ) &
!         call terminate("destroyPETScVars", &
!                        "Could not destroy vector phic")

      call VecDestroy(dJdc, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy vector djdc")
      
      call MatDestroy(dRdWt, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy matrix dRdWt")
      if (ApproxPC) then
         !print *,'mats,prpwpre'
         call MatDestroy(dRdWpret, PETScIerr)
         
         if( PETScIerr/=0 ) &
              call terminate("destroyPETScVars", &
              "Could not destroy matrix dRdWpret")
      endif

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

      call KSPDestroy(ksp, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("destroyPETScVars", &
                       "Could not destroy KSP context")

      ! Synchronize the processors.

      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

#endif

      end subroutine destroyPETScVars
