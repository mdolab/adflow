subroutine saveADjointMatrix(fileName)

  use ADjointPETSc, only: drdwt
  use communication
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"

  ! Input params
  character*(*), intent(in) :: fileName

  ! Working parameters
  PetscViewer binViewer
  integer(kind=intType) :: ierr

  call PetscViewerBinaryOpen(sumb_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatView(dRdwT, binViewer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call PetscViewerDestroy(binViewer,ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine saveADjointMatrix

subroutine saveAdjointPC(fileName)

  use ADjointPETSc, only: drdwpret
  use communication
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"
  ! Input params
  character*(*), intent(in) :: fileName

  ! Working parameters
  PetscViewer binViewer
  integer(kind=intType) :: ierr

  call PetscViewerBinaryOpen(sumb_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatView(dRdwPreT, binViewer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call PetscViewerDestroy(binViewer,ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine saveAdjointPC

subroutine saveADjointRHS(fileName)

  use ADjointPETSc, only: dJdw
  use communication
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"
  ! Input params
  character*(*), intent(in) :: fileName

  ! Working parameters
  PetscViewer binViewer
  integer(kind=intType) :: ierr

  call PetscViewerBinaryOpen(sumb_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecView(dJdw, binViewer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call PetscViewerDestroy(binViewer,ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine saveADjointRHS
