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

      ! Matrices
      call MatDestroy(dRdWT, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      call MatDestroy(dRdWPreT, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      call MatDestroy(dRdx, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      call MatDestroy(dRda, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      call MatDestroy(dFdw, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      call MatDestroy(dFdx, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      ! Vectors

      call VecDestroy(psi, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      call VecDestroy(dJdW, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      call VecDestroy(pvr, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      call VecDestroy(dJdx, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      ! KSP Context
      call KSPDestroy(ksp, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)
#endif

      end subroutine destroyPETScVars
