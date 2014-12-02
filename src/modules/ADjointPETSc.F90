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
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
      ! dRdWT Adjoint matrix that defines the (transpose) linear system dR/dW.
      !       dRdW(il,jl,kl,nw) in each local domain mapped to a matrix.
      !       Size[nNodes*nw, nNodes*nw], where nNodes is the global
      !       number of grid nodes and nw is the number of equations.
      !
      ! psi   Approximate adjoint solution vector solution
      !       of [dRdW]^T {psi} = {dJdw}.
      !       psi(il,jl,kl,nw) in each local domain mapped to a vector.
      !       Size[nNodes*nw].
      !
      ! dJdW  Right-hand side vector holding the cost function Jpetsc
      !       sensitivity w.r.t. W: dJ/dW.
      !       dJdW(il,jl,kl,nw) in each local domain mapped to a vector.
      !       Size[nNodes*nw].
      !            

      ! dRda  Residual R partial sensitivity w.r.t. the design
      !       variable vector alpha: dR/da
      !       dRda(il,jl,kl,nw,n) - d R(i,j,k,nw) / d alpha(n)
      !       in each local domain is mapped to a matrix (m,n)
      !       Size[nNodes*nw,nDesign],where nDesign is the number
      !       of design variables.
      !

      Mat     dRdWT, dRdWPreT, dRdWTShell
      Mat     dRda, dRdx
      Mat     dFcdw, dFcdx, dFndFc
      Mat     dFdx, dFdw, doAdX
      Vec, allocatable, dimension(:,:) :: FMw ! nFM by ntimespectral instance
      Vec, allocatable, dimension(:,:) :: FMx ! nFM by ntimespectral instance
      Mat, allocatable, dimension(:) :: coarsedRdWPreT
      Mat, allocatable, dimension(:) :: restrictionOperator
      Mat, allocatable, dimension(:) :: prolongationOperator
      Vec     psi, dJdW, adjointRHS,adjointRes
      Vec     dJdx, dJdxv
      Vec     gridVec
      Vec     dRdaTPsi, dRdaTPsi_local
      Vec     xVec
      Vec     fVec1, fVec2, overArea, fCell, fNode
      Vec     w_like1, w_like2, psi_like1, psi_like2, x_like
      Vec     dJcdW
      VecScatter dRdaTpsi_scatter
      VecScatter XstoXv
      PetscErrorCode PETScIerr
      PetscFortranAddr   matfreectx(1)

      !adjointKSP   Linear solver (Krylov subspace method) context
      KSP     adjointKSP

      ! KSP Converged Reason
      KSPConvergedReason adjointConvergedReason

      ! Data for dRda and dFMdExtra
      real(kind=realType), allocatable, dimension(:,:) :: dRda_data
      real(kind=realType), allocatable, dimension(:,:,:) :: dFMdExtra
      integer(kind=intType), parameter :: nFM = 8
      integer(kind=intTYpe), parameter:: iSepSensor = 7
      integer(kind=intTYpe), parameter:: iCavitation = 8
      ! Logical identifying the type of PETSc matrix being used for dRdW
       logical :: PETScBlockMatrix
       
      real(kind=realType), allocatable, dimension(:) :: adjResHist
      integer(kind=intType)                          :: adjConvIts
      
      !Binary Viewer
      real(kind=realType) :: localInfo(MAT_INFO_SIZE)
      real(kind=realType) :: sumInfo(MAT_INFO_SIZE)
      real(kind=realType) :: maxInfo(MAT_INFO_SIZE)

#endif
      end module ADjointPETSc
