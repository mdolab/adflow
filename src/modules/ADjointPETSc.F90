!
!     ******************************************************************
!     *                                                                *
!     * File:          ADjointPETSc.F90                                *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 12-15-2005                                      *
!     * Last modified: 06-14-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      module ADjointPETSc
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
#define PETSC_NULL           0
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

      ! dRdW  Adjoint matrix that defines the linear system dR/dW.
      !       dRdW(il,jl,kl,nw) in each local domain mapped to a matrix.
      !       Size[nNodes*nw, nNodes*nw], where nNodes is the global
      !       number of grid nodes and nw is the number of equations.
      !
      ! psi   Approximate adjoint solution vector solution
      !       of [dRdW]^T {psi} = {dJdw}.
      !       psi(il,jl,kl,nw) in each local domain mapped to a vector.
      !       Size[nNodes*nw].
      !
      ! dJdW  Right-hand side vector holding the cost function J
      !       sensitivity w.r.t. W: dJ/dW.
      !       dJdW(il,jl,kl,nw) in each local domain mapped to a vector.
      !       Size[nNodes*nw].
      !            
      ! pvr   Residual vector (r=Ax-b) used to verify the adjoint 
      !       vector solution ( residual = [dRdW]^T {psi} - {dJdw} ).
      !       Size[nNodes*nw].
      ! phic  Structural adjoint vector cast through the mapping onto
      !        the CFD surface. Size [nSurfNodes*3]...Possible changes to 
      !       this defintion in the future
	

      Mat     dRdW ,dRdWFD,dRdWPreCond
      Vec     psi, dJdW, pvr,phic,dJcdW

      ! ksp   Linear solver (Krylov subspace method) context
      ! pc    Preconditioner context

      KSP     ksp
      PC      pc

      !Subcontexts for ASM preconditioner
      PC subpc
      !KSP, dimension(:),allocatable :: subksp!(*)! = PETSC_NULL
      KSP subksp

      !index sets for the ASM preconditioner
      !IS,    dimension(10) ::  is
      IS  is

      !Additional variables for the ASM preconditioner
      PetscInt nlocal,first,Nsub,length
      !PetscInt ,    dimension(10) :: idx

      ! dRda  Residual R partial sensitivity w.r.t. the design
      !       variable vector alpha: dR/da
      !       dRda(il,jl,kl,nw,n) - d R(i,j,k,nw) / d alpha(n)
      !       in each local domain is mapped to a matrix (m,n)
      !       Size[nNodes*nw,nDesign],where nDesign is the number
      !       of design variables.

      Mat     dRda, dRdsigma, dRdx,dRdxFD, dRdy, dRdz,dSdw,dSdx,dSdxfd2
      Mat     dCda, dCdx, dCdw
      Vec     dJdC
!###      ! Variables for the insertion of dRdx
!###      ! row           Row number
!###      ! row_size      Number of elements in each row of dRdx
!###      ! column(n)     Column indices for each row
!###      ! row_values(n) Values for each column of the given row
!###
!###      PetscInt    row, row_size, column(100)
!###      PetscReal   row_values(100)
!###
      ! dJda  Cost function J partial sensitivity w.r.t. the design
      !       variable vector alpha: dJ/da. Row vector.
      !       Size[nDesign].
      !
      ! dIda  Cost function I total sensitivity w.r.t. the design
      !       variable vector alpha: dI/da, computed as
      !       {dIdW} = {dJdW} - {psi}^T [dRda]. Row vector.
      !       Size[nDesign].

      Vec     dJda, dIda 
      Vec     dJdx, dJcdx,dJdy, dJdz, dIdx, dIdy, dIdz,dIdxs

      ! Some useful real constants

      PetscScalar   PETScNegOne, PETScZero, PETScOne

      ! Some auxiliar integer variables

      PetscInt      PETScIntM, PETScIntN

      ! Some auxiliar integer variables

      PetscScalar   PETScScalarV

      ! Some variables for debugging

      PetscInt columnsize
      PetscInt columnindex(625)
      PetscScalar columnvalues(625)
      PetscInt row, column
      PetscScalar value
      PetscScalar PetscSubt
      PetscScalar PetscAdd	
      PetscScalar ILULevels 

      ! Logical identifying the type of PETSc matrix being used for dRdW

      logical :: PETScBlockMatrix

!moved to inputADjoint
!!$      ! Convergence tolerances and monitor parameters used by PETSc to
!!$      ! solve the adjoint equations. These are set by the routine
!!$      ! "KSPSetTolerances" and "KSPGMRESSetRestart".
!!$      !
!!$      ! adjRelTol    Relative tolerance
!!$      ! adjAbsTol    Absolute tolerance
!!$      ! adjDivTol    Relative tolerance increase to divergence
!!$      ! adjMaxIter   Maximum number of iterations
!!$      ! adjRestart   Maximum number of steps before restart
!!$      !              It has a high impact on the required memory!
!!$      ! adjMonStep   Convergence monitor step
!!$
!!$      real(kind=realType),   parameter :: adjRelTol  = 1.0e-12_realType
!!$      real(kind=realType),   parameter :: adjAbsTol  = 1.0e-25_realType
!!$      real(kind=realType),   parameter :: adjDivTol  = 1.0e+5_realType
!!$      integer(kind=intType), parameter :: adjMaxIter = 500_intType!128_intType
!!$      integer(kind=intType), parameter :: adjRestart = 80_intType!128_intType
!!$      integer(kind=intType), parameter :: adjMonStep = 10_intType

      ! adjResHist   Linear solver residual history [adjMaxIter]
      ! adjConvIts   Number of iterations to convergence

      real(kind=realType), allocatable, dimension(:) :: adjResHist
      integer(kind=intType)                          :: adjConvIts

      !scalars for addition

      !Varibles for accessing the array psi and psidir
       
      PetscScalar x_array(1)

      PetscOffset i_x

#endif

      end module ADjointPETSc
