!
!     ******************************************************************
!     *                                                                *
!     * This file contains three routines for saving the adjoint       *
!     * matrix, the preconditioning matrix and the RHS to a file.      *
!     * The intended usage for this is to do experimentation with      *
!     * stand-alone PETSc programs that just solve a linear system.    *
!     *                                                                *
!     ******************************************************************

subroutine saveADjointMatrix(fileName)

  use constants
  use ADjointPETSc, only: drdwt
  use communication, only : sumb_comm_world
  use utils, only : EChk
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

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

  use constants
  use ADjointPETSc, only: drdwpret
  use communication, only : sumb_comm_world
  use utils, only : EChk
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

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

!   use ADjointPETSc, only: dJdw
!   use communication
!   implicit none

! #ifndef USE_NO_PETSC
! #define PETSC_AVOID_MPIF_H
! #include "finclude/petsc.h"
!   ! Input params
!   character*(*), intent(in) :: fileName

!   ! Working parameters
!   PetscViewer binViewer
!   integer(kind=intType) :: ierr

!   call PetscViewerBinaryOpen(sumb_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   call VecView(dJdw, binViewer, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   call PetscViewerDestroy(binViewer,ierr)
!   call EChk(ierr, __FILE__, __LINE__)
! #endif
end subroutine saveADjointRHS

subroutine saveCellCenters(fileName)

  use blockPointers
  use iteration
  use inputTimeSpectral
  use adjointVars, only: nCellsLocal
  use communication
  use utils, only : setPointers, EChk
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  ! Input params
  character*(*), intent(in) :: fileName

  ! Working parameters
  PetscViewer binViewer
  integer(kind=intType) :: nn,sps,i,j,k,n
  integer(kind=intType) :: ierr,nDimC, iRow,level
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  real(kind=realType),dimension(3)::cellCenter
  integer(kind=intType),dimension(3)::iCol
  
  Mat cellCenters

  nDimC = nCellsLocal(1_intType)*nTimeIntervalsSpectral
  ! create a temporary matrix for the cell centers
  allocate(nnzDiagonal(nDimC), nnzOffDiag(nDimC))
  do i=1,nDimC
     nnzDiagonal(i)=3
     nnzOffDiag(i)=3
  end do

  call myMatCreate(cellCenters, 1, nDimC, 3, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)
  deallocate(nnzDiagonal, nnzOffDiag)
  
  iCol(1) = 0
  iCol(2) = 1
  iCOl(3) = 2
  ! compute and store the cell centers
  level = currentLevel
  domainLoop: do nn=1, nDom
     spectralLoop: do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        do k=1,kl
           do j=1,jl
              do i=1,il
                 iRow = flowDoms(nn, level, sps)%globalCell(i, j, k)
                 ! The location of the cell center is determined
                 ! by averaging the cell coordinates.
                 do n=1,3
                    cellCenter(n) = (x(i-1,j-1,k-1,n) + x(i,j-1,k-1,n)  &
                         +  x(i-1,j,  k-1,n) + x(i,j,  k-1,n)  &
                         +  x(i-1,j-1,k,  n) + x(i,j-1,k,  n)  &
                         +  x(i-1,j,  k,  n) + x(i,j,  k,  n))/8
                 end do
                 call MatSetValues(cellCenters, 1, iRow, 3, iCol, cellCenter, &
                      INSERT_VALUES, ierr)
              end do
           end do
        end do
     end do spectralLoop
  end do domainLoop
  
! PETSc Matrix Assembly begin
  call MatAssemblyBegin(cellCenters, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  ! Complete the matrix assembly.
  call MatAssemblyEnd  (cellCenters, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call PetscViewerBinaryOpen(sumb_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatView(cellCenters, binViewer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call PetscViewerDestroy(binViewer,ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine saveCellCenters
