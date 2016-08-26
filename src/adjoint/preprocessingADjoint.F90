subroutine preprocessingADjoint
  !
  !      Perform the preprocessing tasks for the adjoint solver. This   
  !      routine is called only once. The memory allcoated here is      
  !      deallocated in src/utils/releaseMemory.f90                     
  !
  use constants
  use communication, only : sumb_comm_world
  use adjointVars, only :nCellsLocal, nNOdesLocal
  use flowVarRefState, only : nw, nwf
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use inputAdjoint, only : frozenTurbulence
  use ADjointPETSc, only: w_like1, w_like2, PETScIerr, &
       psi_like1, psi_like2, x_like, psi_like3
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

#endif

  !     Local variables.
  !
  integer(kind=intType) :: ndimW, ncell, nState, nDimPsi, nDimX
  !
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

  nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
  nDimPsi = nState*  nCellsLocal(1_intType)*nTimeIntervalsSpectral
  nDimX = 3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral

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

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,nstate,ndimPsi,PETSC_DECIDE, &
       PETSC_NULL_SCALAR,psi_like2,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,nstate,ndimPsi,PETSC_DECIDE, &
       PETSC_NULL_SCALAR,psi_like3,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,3,ndimX,PETSC_DECIDE, &
       PETSC_NULL_SCALAR,x_like,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! Need to initialize the stencils as well, only once:
  call initialize_stencils

end subroutine preprocessingADjoint
