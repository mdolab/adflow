subroutine setupAdjointMatrix
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Assembles the adjoint (dRdw)^T Matrix. Called from Python      *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use inputADjoint
  use communication
  implicit none

  ! Local variables.
  logical :: useAD, useTranspose, usePC, useObjective
  integer(kind=intType) :: ierr
  real(kind=realType), dimension(2) :: time
  real(kind=realType)               :: timeAdjLocal, timeAdj
  interface
     subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
          useObjective, frozenTurb, level, matrixTurb)
       use precision
       implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
       Mat :: matrix
       Mat, optional :: matrixTurb
       ! Input Variables
       logical, intent(in) :: useAD, usePC, useTranspose, useObjective, frozenTurb
       integer(kind=intType), intent(in) :: level
     end subroutine setupStateResidualMatrix
  end interface

#ifndef USE_NO_PETSC

  if( myid ==0 ) &
       write(*, 10) "Assembling State Residual Matrices..."

  call cpu_time(time(1))

  ! Compute dRdw with forward AD
  useAD = .True.
  usePC = .False.
  useTranspose = .True.
  useObjective = .True.

  call setupStateResidualMatrix(drdwt, useAD, usePC, useTranspose, useObjective, &
       frozenTurbulence, 1_intType)

  call cpu_time(time(2))
  timeAdjLocal = time(2)-time(1)

  ! Reudce and display max
  call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
       mpi_max, 0, SUMB_COMM_WORLD, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)

  if(myid ==0)  then
     write(*, 20) "Assembling State Residual Matrices Time (s) = ", timeAdj
  end if

  ! Output formats.
10 format(a)
20 format(a, 1x, f8.2)

#endif
end subroutine setupAdjointMatrix
