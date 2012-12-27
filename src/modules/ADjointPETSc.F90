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
#include "include/petscversion.h"

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
      ! dJdW  Right-hand side vector holding the cost function J
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

      Mat     dRdWT,dRdWPreT, matTemp, matTemp2
      Mat     dRda, dRdx,dFdw,dFdx
      Vec     psi, dJdW, adjointRHS,adjointRes
      Vec     dJdx
      Vec     gridVec
      Vec     dRdaTPsi, dRdaTPsi_local
      Vec     xVec
      Vec     fVec1, fVec2
      Vec     w_like1, w_like2
      Vec     dJcdW
      ! ksp   Linear solver (Krylov subspace method) context
      ! pc    Preconditioner context
      VecScatter dRdaTpsi_scatter
      
      KSP     ksp, master_PC_KSP
      PC      pc, master_PC

      !Subcontexts for ASM preconditioner
      PC subpc
      KSP subksp

      ! KSP Converged Reason
      KSPConvergedReason adjointConvergedReason

      !index sets for the ASM preconditioner
      IS  is

      PetscFortranAddr,dimension(:),allocatable :: is_local
      PetscFortranAddr,dimension(:),allocatable :: sub_ksp_list

      !Additional variables for the ASM preconditioner
      PetscInt nlocal,first,Nsub,length
   
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

      ! Data for dRda
      real(kind=realType), allocatable, dimension(:,:) :: dRda_data

      ! Logical identifying the type of PETSc matrix being used for dRdW

      logical :: PETScBlockMatrix

      real(kind=realType), allocatable, dimension(:) :: adjResHist
      integer(kind=intType)                          :: adjConvIts

      !scalars for addition

      !Varibles for accessing the array psi and psidir
       
      PetscScalar x_array(1)

      PetscOffset i_x

      !Binary Viewer
      PetscViewer Bin_Viewer
      real(kind=realType) :: localInfo(MAT_INFO_SIZE)
      real(kind=realType) :: sumInfo(MAT_INFO_SIZE)
      real(kind=realType) :: maxInfo(MAT_INFO_SIZE)

#endif

      end module ADjointPETSc
