!
!     ******************************************************************
!     *                                                                *
!     * File:          getADjoint.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 08-18-2008                                      *
!     * Last modified: 10-04-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine getADjoint(ncells,functionGradLocal)
  use ADjointPETSc

  implicit none
  !
  !     Subroutine arguments.
  !
  integer(kind=intType),intent(in):: ncells
  real(kind=realType),dimension(ncells),intent(out) :: functionGradLocal

  ! Local Variables
  integer(kind=intType) :: idxmg, iLow, iHigh, n 		

#ifndef USE_NO_PETSC		

  !     ******************************************************************
  !     *                                                                *
  !     * Transfer solution from PETSc context.                          *
  !     *                                                                *
  !     ******************************************************************
  !

  ! Query about the ownership range.
  ! iHigh is one more than the last element stored locally.

  call VecGetOwnershipRange(psi, iLow, iHigh, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  n = 0
  do idxmg=iLow, iHigh-1
     n = n + 1
     call VecGetValues(psi, 1, idxmg, &
          functionGradLocal(n), PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
  enddo

#endif

end subroutine getADjoint
