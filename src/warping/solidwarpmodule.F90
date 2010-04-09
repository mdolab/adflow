      module solidwarpmodule
!
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

      PetscErrorCode PETScIerr
      PetscMPIInt    PETScRank, PETScSize
      
      Vec     us
      Vec     uu
      Vec     fu
      Mat     Kuu,Kus
      IS      row_perm,col_perm
      KSP     ksp
      PC      pc
      PetscScalar   PETScNegOne, PETScZero, PETScOne
      MatFactorInfo  info(MAT_FACTORINFO_SIZE)

      integer(kind=intType),dimension(:), allocatable :: &
           allNDom,cumNDom,allNElem,cumNElem,allBlockIDs
      integer(kind=intType),dimension(:),allocatable :: &
           Kdispls,BCdispls,Krecvcount,BCrecvcount
      integer(kind=intType) :: nElem
      
      ! Storage for stiffness matrices/Boundary Conditions
      real(kind=realType),  dimension(:), allocatable :: &
           allK,localK,allBCVal,localBCVal
      real(kind=realType),  dimension(:,:,:,:), allocatable:: Xtarget


#endif

    end module solidwarpmodule
