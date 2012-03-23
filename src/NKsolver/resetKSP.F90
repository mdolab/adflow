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

!   if (NKSolverSetup) then

!      call KSPReset(global_ksp)

!   end if
end subroutine resetKSP


      
