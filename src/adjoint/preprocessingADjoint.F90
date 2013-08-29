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
  !     * Perform the preprocessing tasks for the adjoint solver         *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use adjointVars
  use flowVarRefState
  use inputTimeSpectral
  use inputAdjoint
  use ADjointPETSc, only: w_like1, w_like2, fvec1, fvec2, PETScIerr, &
       psi_like1, psi_like2
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#endif

  !     Local variables.
  !
  integer(kind=intType) :: ierr, level, ndimW, ndimS, nTS, ncell, nState, nDimPsi
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
  ! Create PETSc Vectors that are actually empty. These do NOT take
  ! any (substantial) memory. We want to keep these around inbetween
  ! creations/deletions of adjoint/NKsolver memory

  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  endif

  ! Create two (empty) Vectors for getdFdx(T)Vec operations
  call getForceSize(nDimS,ncell)
  nDimS = nDimS * 3 *nTimeIntervalsSpectral! Multiply by 3 for each
                                           ! dof on each point

  nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
  nDimPsi = nState*  nCellsLocal(1_intType)*nTimeIntervalsSpectral

  if (.not. adjointInitialized) then
     if (PETSC_VERSION_MINOR == 2) then
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndimS,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,fVec1,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndimS,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,fVec2,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        ! Two w-like vectors. 
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndimW,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,w_like1,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndimW,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,w_like2,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)

        ! Two psi-like vectors. 
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndimPsi,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,psi_like1,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndimPsi,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,psi_like2,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)

     else
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,1,ndimS,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,fVec1,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,1,ndimS,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,fVec2,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        ! Two w-like vectors. 
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,nw,ndimW,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,w_like1,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,nw,ndimW,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,w_like2,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        ! Two psi-like vectors. 
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,nState,ndimPsi,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,psi_like1,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,nState,ndimPsi,PETSC_DECIDE, &
             PETSC_NULL_SCALAR,psi_like2,PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)


     end if
     adjointInitialized = .True.
  endif
  ! Need to initialize the stencils as well, only once:
  call initialize_stencils
     
end subroutine preprocessingADjoint
