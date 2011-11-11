!
!     ******************************************************************
!     *                                                                *
!     * File:          setADjoint.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 10-04-2010                                      *
!     * Last modified: 10-04-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine setADjoint(ncells,functionGradLocal)
  use ADjointPETSc

  implicit none
  !
  !     Subroutine arguments.
  ! 
  integer(kind=intType),intent(in):: ncells
  real(kind=realType),dimension(ncells),intent(in) :: functionGradLocal
  !
  !     Local variables.
  integer(kind=intType) :: idxmg, iLow, iHigh, n 		

#ifndef USE_NO_PETSC		
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Transfer stored adjoint to PETSc                               *
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
     call VecSetValue(psi, idxmg, &
          functionGradLocal(n),INSERT_VALUES, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
  enddo

  call VecAssemblyBegin(psi, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  call VecAssemblyEnd(psi,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

#endif
end subroutine setADjoint
