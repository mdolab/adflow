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
subroutine getADjoint(nstate,adjoint)

#ifndef USE_NO_PETSC	
#define PETSC_AVOID_MPIF_
#include "finclude/petscdef.h"

  use ADjointPETSc, only : psi
  use petscvec
  use constants

  implicit none
  !
  !     Subroutine arguments.
  !
  integer(kind=intType),intent(in):: nstate
  real(kind=realType),dimension(nstate),intent(out) :: adjoint
  real(kind=realType),pointer :: psi_pointer(:)

  ! Local Variables
  integer(kind=intType) :: i, ierr

  ! Copy out adjoint vector:
  call VecGetArrayF90(psi,psi_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Do a straight copy:
  do i=1,nstate
     adjoint(i) = psi_pointer(i)
  end do

  call VecRestoreArrayF90(psi,psi_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

end subroutine getADjoint

