!
!     ******************************************************************
!     *                                                                *
!     * File:          createPETScKsp.F90                              *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader                *
!     * Starting date: 12-15-2005                                      *
!     * Last modified: 01-19-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine createPETScKsp
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the Krylov subspace linear solver context,              *
  !     * the preconditioner context and set their various options.      *
  !     * This defines the linear solver to be used to solve the adjoint *
  !     * system of equations.                                           *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  implicit none

  
#ifndef USE_NO_PETSC

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the ksp context.                                        *
  !     *                                                                *
  !     ******************************************************************

  call KSPCreate(SUMB_PETSC_COMM_WORLD, ksp, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#endif

end subroutine createPETScKsp

