module solidwarpmodule
  
  use constants
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H
  !
  !     ******************************************************************
  !     *                                                                *
  !     * PETSc include files.                                           *
  !     *                                                                *
  !     * This module uses CPP for preprocessing, as indicated by the use*
  !     * of PETSc include files in the directory petsc/include/finclude.*
  !     * This convention enables use of the CPP preprocessor, which     *
  !     * allows the use of the #include statements that define PETSc    *
  !     * objects and variables.                                         *
  !     *                                                                *
  !     * Include statements are required for KSP Fortran programs:      *
  !     *                                                                *
  !     *   petsc.h        Base PETSc routines                           *
  !     *   petscvec.h     Vectors                                       *
  !     *   petscvec.h90   Vectors (additional F90 pointers support)     *
  !     *   petscmat.h     Matrices                                      *
  !     *   petscmat.h90   Matrices (additional F90 pointers support)    *
  !     *   petscksp.h     Krylov subspace methods                       *
  !     *   petscpc.h      Preconditioners                               *
  !     *   petscviewer.h  Viewers                                       *
  !     *   petscis.h      Index sets                                    *
  !     *   petscis.h90    Index sets (additional F90 pointers support)  *
  !     *   petscdraw.h    PETSc draw                                    *
  !     *   petscmg.h      MG preconditioner                             *
  !     *   petscsys.h     System package                                *
  !     *                                                                *
  !     ******************************************************************
  !
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscsys.h"

  ! -----------------
  ! PETSC Variables
  ! -----------------
  Vec     us    ! structural state vector,u, specified part, size nus
  Vec     uu    ! structural state vector,u, unknown part, size nuu
  Vec     u     ! full structural vector u=[uu us]^T
  Vec     fu    ! structural rhs force vector, size nuu
  Vec     dXv   ! A vector containing the current derivative of
  ! volume coords wrt surface coord

  ! K matrix is partitioned according to:
  ! K = [kuu kus]
  !     [ksu kss]
  Mat     Kuu,Kus !kuu is unknown-unknown part and Kus is unknown-known part
  Mat     dXvdXvFE ! derivative of the volume coordinate wrt all
  ! structural u. This is just the interpolation shape functions


  !Vec     us_local ! A (full) local copy of us
  !Vec     uu_local ! A (full) local copy of uu
  Vec      u_local ! A (full) local copy of u 
  VecScatter u_scatter
  !VecScatter uu_scatter,us_scatter ! Scatter contexts for
  ! scattering uu and us to all
  ! procs
!  IS      row_perm,col_perm

  KSP     ksp
  PC      pc

  ! Define some common petsc scalars
  PetscScalar   PETScNegOne, PETScZero, PETScOne
  MatFactorInfo  info(MAT_FACTORINFO_SIZE)

  ! -----------------------------
  ! Variables set from Python
  ! ----------------------------
  integer(kind=intType),dimension(:),allocatable :: &
       l_index,lptr,md_g_index,md_g_ptr
   integer(kind=intType),dimension(:,:),allocatable :: &
        l_sizes

  integer(kind=intType) :: nuu,nus,nblock,nli,nuu_local,nus_local

  integer(kind=intType),dimension(:), allocatable :: uu_ownership,uu_indices
  integer(kind=intType),dimension(:), allocatable :: us_ownership,us_indices

 !  ! ----------------------------------------------
!   ! Variables for Serial Version of Mesh Warping
!   ! ----------------------------------------------
!   integer(kind=intType),dimension(:), allocatable :: &
!        allNDom,cumNDom,allNElem,cumNElem,allBlockIDs
!   integer(kind=intType),dimension(:),allocatable :: &
!        Kdispls,BCdispls,Krecvcount,BCrecvcount
!   real(kind=realTYpe),dimension(:),allocatable :: &
!        allK,allBCVal
!   integer(kind=intType) :: nElem

#endif

end module solidwarpmodule
