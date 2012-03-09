!
!     ******************************************************************
!     *                                                                *
!     * File:          resetKSP.F90                                    *
!     * Author:        Gaetan Kenway                                   *
!     *                                                                *
!     ******************************************************************
!
subroutine resetKSP

  ! Reset the KSP object used for the NK solver
       
  use NKSolverVars
  use communication 
  implicit none

  ! Local Vars
  integer(kind=intType) :: ierr, nlocal, first
  logical :: useAD, usePC, useTranspose

  if (NKSolverSetup) then

     ! Destroy Solver
     call KSPDestroy(global_ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Create
     call KSPCreate(SUMB_PETSC_COMM_WORLD,global_ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set operators for the solver
     call KSPSetOperators(global_ksp,dRdw,dRdWPre, &
          DIFFERENT_NONZERO_PATTERN,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call setupNKKSP()

  end if
end subroutine resetKSP


      
