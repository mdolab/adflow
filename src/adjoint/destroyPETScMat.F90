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
     
      call MatDestroy(dRdWT, PETScIerr)
      call EChek(PETScIerr,__file__,__line__)

      call MatDestroy(dRdWPreT, PETScIerr)
      call EChek(PETScIerr,__file__,__line__)

      call MatDestroy(dRdx, PETScIerr)
      call EChek(PETScIerr,__file__,__line__)

      call MatDestroy(dRda, PETScIerr)
      call EChek(PETScIerr,__file__,__line__)

      call MatDestroy(dFdw, PETScIerr)
      call EChek(PETScIerr,__file__,__line__)

      call MatDestroy(dFdx, PETScIerr)
      call EChek(PETScIerr,__file__,__line__)
 
#endif

    end subroutine destroyPETScMat
