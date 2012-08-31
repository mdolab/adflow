!
!     ******************************************************************
!     *                                                                *
!     * File:          preprocessingADjoint.f90                        *
!     * Author:        Andre C. Marta , C.A.(Sandy) Mader              *
!     *                Seongim Choi
!     * Starting date: 07-20-2006                                      *
!     * Last modified: 01-18-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine preprocessingADjoint
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Perform the preprocessing tasks for the adjoint solver, such   *
  !     * as assertions testing, memory allocation and global indexing.  *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use adjointVars
  use precision
  use flowVarRefState
  use inputTimeSpectral
  use ADjointPETSc, only: w_like1, w_like2, fvec1, fvec2, PETScIerr
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#endif

  !     Local variables.
  !
  integer(kind=intType) :: ierr,level,ndimW,ndimS,nTS
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !

  ! Test the discrete adjoint solver assertions.
  level = 1_intType

  call assertionsADjoint(level)
  
  ! Allocate the memory for the global cell indexing of the
  ! computational mesh, later used to assemble the global
  ! adjoint system of equations.

  call allocMemADjoint(level)
     
  ! Determine the global cell and Node numbering.
  call setGlobalCellsAndNodes(level)

  ! Create PETSc Vectors that are actually empty. These do NOT take
  ! any (substantial) memory. We want to keep these around inbetween
  ! creations/deletions of adjoint/NKsolver memory

  ! Create two (empty) Vectors for getdFdx(T)Vec operations
  call getForceSize(nDimS,nTS)
  nDimS = nDimS * 3 *nTimeIntervalsSpectral! Multiply by 3 for each
                                           ! dof on each point

  nDimW = nw * nCellsLocal*nTimeIntervalsSpectral

  if (PETSC_VERSION_MINOR == 2) then
     call VecCreateMPIWithArray(SUMB_PETSC_COMM_WORLD,ndimS,PETSC_DECIDE, &
          PETSC_NULL_SCALAR,fVec1,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)

     call VecCreateMPIWithArray(SUMB_PETSC_COMM_WORLD,ndimS,PETSC_DECIDE, &
          PETSC_NULL_SCALAR,fVec2,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)

     ! Two w-like vectors. 
     call VecCreateMPIWithArray(SUMB_PETSC_COMM_WORLD,ndimW,PETSC_DECIDE, &
          PETSC_NULL_SCALAR,w_like1,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call VecCreateMPIWithArray(SUMB_PETSC_COMM_WORLD,ndimW,PETSC_DECIDE, &
          PETSC_NULL_SCALAR,w_like2,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
  else
     call VecCreateMPIWithArray(SUMB_PETSC_COMM_WORLD,1,ndimS,PETSC_DECIDE, &
          PETSC_NULL_SCALAR,fVec1,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)

     call VecCreateMPIWithArray(SUMB_PETSC_COMM_WORLD,1,ndimS,PETSC_DECIDE, &
          PETSC_NULL_SCALAR,fVec2,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)

     ! Two w-like vectors. 
     call VecCreateMPIWithArray(SUMB_PETSC_COMM_WORLD,nw,ndimW,PETSC_DECIDE, &
          PETSC_NULL_SCALAR,w_like1,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call VecCreateMPIWithArray(SUMB_PETSC_COMM_WORLD,nw,ndimW,PETSC_DECIDE, &
          PETSC_NULL_SCALAR,w_like2,PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
  end if
     
end subroutine preprocessingADjoint
