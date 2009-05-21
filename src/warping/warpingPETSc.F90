!
!     ******************************************************************
!     *                                                                *
!     * File:          warpingPETSc.F90                                *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 05-21-2009                                      *
!     * Last modified: 05-21-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      module warpingPETSc
!
!     ******************************************************************
!     *                                                                *
!     * This module contains the objects used by PETSc for the         *
!     * solution of the discrete adjoint equations.                    *
!     *                                                                *
!     ******************************************************************
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
!#f90##include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
!#f90##include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
!#include "include/finclude/petscviewer.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
!#include "include/finclude/petscdraw.h"
!#include "include/finclude/petscmg.h"
!#include "include/finclude/petscsys.h"

!
!     ******************************************************************
!     *                                                                *
!     * PETSc variables declaration.                                   *
!     *                                                                *
!     ******************************************************************
!
     
      ! PETScIerr  Error return code
      ! PETScRank  Processor number (0-based index)
      ! PETScSize  Number of processors

      PetscErrorCode PETScIerr
      PetscMPIInt    PETScRank, PETScSize

      ! dXvdXs Matrix of derivatives of the volume mesh coordinates with
      !        respect to the surface coordinates.
      !       Size[nNodes*3, nNodesSurf*3], where nNodes is the global
      !       number of grid nodes and nNodesSurf is the number surface
      !       nodes.
      !
 
      Mat     dXvdXs ,dXvdXsFD

#endif

      end module warpingPETSc
