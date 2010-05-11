!
!     ******************************************************************
!     *                                                                *
!     * File:          createPETScKsp.F90                              *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader                *
!     * Starting date: 12-15-2005                                      *
!     * Last modified: 01-19-2010                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine createPETScKsp(level)
!
!     ******************************************************************
!     *                                                                *
!     * Create the Krylov subspace linear solver context,              *
!     * the preconditioner context and set their various options.      *
!     * This defines the linear solver to be used to solve the adjoint *
!     * system of equations.                                           *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars
      use blockPointers
      use flowVarRefState  !nw
      use inputADjoint
      implicit none

!!$!
!!$!     User-defined Fortran routine.
!!$!
!!$      external MyKSPMonitor

!
!     Subroutine Variables
!

      integer(kind=intType)::level

!
!     Local variables.
!
      real(kind=realType)   :: rTol, aTol, dTol
      integer(kind=intType) :: mIts
      character(len=10)     :: kspType, pcType

#ifndef USE_NO_PETSC

!!$      ! PETSc macros are lost and have to be redefined.
!!$      ! They were extracted from: <PETSc root dir>/include/petscksp.h
!!$      !                                                   /petscpc.h
!!$!Solvers
!!$#define KSPGMRES      "gmres"   
!!$#define KSPBCGS       "bcgs"
!!$#define KSPCG         "cg"
!!$#define KSPFGMRES     "fgmres"
!!$
!!$!Global Preconditioners
!!$#define PCJACOBI      "jacobi"
!!$#define PCBJACOBI     "bjacobi"
!!$#define PCASM         "asm"
!!$
!!$!Local Preconditioners
!!$#define PCILU         "ilu"
!!$#define PCICC         "icc"
!!$#define PCLU          "lu"
!!$#define PCCHOLESKY    "cholesky"
!!$
!!$!Matrix Reorderings
!!$#define MATORDERING_NATURAL       "natural"
!!$#define MATORDERING_RCM       "rcm"
!!$#define MATORDERING_ND        "nd"
!!$#define MATORDERING_1WD       "1wd"
!!$#define MATORDERING_QMD       "qmd"
!!$
!!$!Other miscellaneous definitions
!!$#define PETSC_NULL           0
!!$#define PETSC_DEFAULT        -2
!!$#define KSPPREONLY    "preonly"


!
!     ******************************************************************
!     *                                                                *
!     * Create the ksp context.                                        *
!     *                                                                *
!     ******************************************************************
!
      ! KSPCreate - Creates the default KSP context.
      !
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPCreate(MPI_Comm comm,KSP *inksp, PetscErrorCode ierr)
      !
      ! Collective on MPI_Comm
      !
      ! Input Parameter
      !   comm - MPI communicator
      !
      ! Output Parameter
      !   ksp  - location to put the KSP context
      !
      ! Notes
      ! The default KSP type is GMRES with a restart of 30, using
      !   modified Gram-Schmidt orthogonalization.
      !
      ! see .../petsc/docs/manualpages/KSP/KSPCreate.html
      ! or PETSc users manual, pp.63

      call KSPCreate(PETSC_COMM_WORLD, ksp, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPCreate")  
!!$!
!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     * Set matrix to provide the values for the KSP operator          *
!!$!     *                                                                *
!!$!     ******************************************************************
!!$!
!!$      ! KSPSetOperators - Sets the matrix associated with the linear
!!$      !                   system and a (possibly) different one
!!$      !                   associated with the preconditioner.
!!$      ! Synopsis
!!$      !
!!$      ! #include "petscksp.h" 
!!$      ! call KSPSetOperators(KSP ksp,Mat Amat,Mat Pmat, &
!!$      !                      MatStructure flag,PetscErrorCode ierr)
!!$      !
!!$      ! Collective on KSP and Mat
!!$      !
!!$      ! Input Parameters
!!$      !   ksp  - the KSP context
!!$      !   Amat - the matrix associated with the linear system
!!$      !   Pmat - the matrix to be used in constructing the 
!!$      !          preconditioner, usually the same as Amat.
!!$      !   flag - flag indicating information about the preconditioner
!!$      !          matrix structure during successive linear solves.
!!$      !          This flag is ignored the first time a linear system is
!!$      !          solved, and thus is irrelevant when solving just one
!!$      !          linear system.
!!$      !
!!$      ! see .../petsc/docs/manualpages/KSP/KSPSetOperators.html
!!$      ! or PETSc users manual, pp.63
!!$
!!$      if (ApproxPC)then
!!$         
!!$         !setup the approximate PC Matrix
!!$         !call setupADjointPCMatrix(level)
!!$         call setupADjointPCMatrixTranspose(level)
!!$         
!!$         !now set up KSP Context
!!$         !call KSPSetOperators(ksp,dRdW,dRdWPre, &
!!$         !                  DIFFERENT_NONZERO_PATTERN,PETScIerr)
!!$         call KSPSetOperators(ksp,dRdWT,dRdWPreT, &
!!$              DIFFERENT_NONZERO_PATTERN,PETScIerr)
!!$
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKSP", "Error in KSPSetOperators.")
!!$
!!$
!!$
!!$         !call PetscOptionsPrint(PETScIerr)
!!$         !call PetscOptionsSetValue('-ksp_gmres_modifiedgramschmidt',PETSC_NULL_CHARACTER,PETScIerr)
!!$         
!!$         !call PetscOptionsPrint(PETScIerr)
!!$
!!$
!!$      else
!!$         
!!$         ! Use the exact jacobian.
!!$         ! Here the matrix that defines the linear system
!!$         ! also serves as the preconditioning matrix.
!!$         
!!$         !call KSPSetOperators(ksp,dRdW,dRdW, &
!!$         !     DIFFERENT_NONZERO_PATTERN,PETScIerr)
!!$         call KSPSetOperators(ksp,dRdWT,dRdWT, &
!!$              DIFFERENT_NONZERO_PATTERN,PETScIerr)
!!$         
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKSP", "Error in KSPSetOperators.")
!!$
!!$      end if
!!$
!!$!
!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     * Set the various options for the KSP.                           *
!!$!     *                                                                *
!!$!     ******************************************************************
!!$!
!!$      ! KSPSetFromOptions - Sets KSP options from the options database.
!!$      !   This routine must be called before KSPSetUp() if the user is
!!$      !   to be allowed to set the Krylov type.
!!$      !
!!$      ! Synopsis
!!$      !
!!$      ! #include "petscksp.h" 
!!$      ! call KSPSetFromOptions(KSP ksp, PetscErrorCode ierr)
!!$      !
!!$      ! Collective on KSP
!!$      !
!!$      ! Input Parameters
!!$      !   ksp - the Krylov space context
!!$      !
!!$      ! see .../petsc/docs/manualpages/KSP/KSPSetFromOptions.html
!!$
!!$      call KSPSetFromOptions(ksp, PETScIerr)
!!$
!!$      if( PETScIerr/=0 ) &
!!$        call terminate("createPETScKsp", "Error in KSPSetFromOptions")
!!$
!!$      !call KSPSetComputeSingularValues(ksp, PETSC_TRUE, PETScIerr)
!!$
!!$
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Now setup the specific options for the selected solver        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$      
!!$      select case(ADjointSolverType)
!!$         
!!$      case(PETSCGMRES)
!!$      
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set the Krylov method.                                        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$         
!!$         ! KSPSetType - Builds KSP for a particular solver.
!!$         ! 
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPSetType(KSP ksp, KSPType type, PetscErrorCode ierr)
!!$         !
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameters
!!$         !   ksp  - the Krylov space context
!!$         !   type - a known method: KSPGMRES - GMRES (default)
!!$         !                          KSPCG    - Conjugate gradients
!!$         !                          KSPLSQR  -Least-squares method
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPSetType.html
!!$         ! or PETSc users manual, pp.65
!!$         
!!$         call KSPSetType(ksp, KSPGMRES, PETScIerr)
!!$         
!!$         
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", "Error in KSPSetType")
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Get the Krylov method.                                        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$         ! KSPGetType - Gets the KSP type as a string from the KSP object.
!!$         ! 
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPGetType(KSP ksp, KSPType type, PetscErrorCode ierr)
!!$         !
!!$         ! Not collective
!!$         !
!!$         ! Input Parameters
!!$         !   ksp  - the Krylov space context
!!$         ! Output Parameters
!!$         !   type - KSP type
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPGetType.html
!!$         
!!$         if( PETScRank==0 ) then
!!$            call KSPGetType(ksp, kspType, PETScIerr)
!!$            if( PETScIerr/=0 ) &
!!$                 call terminate("createPETScKsp", "Error in KSPGetType")
!!$            write(*,10) kspType
!!$         endif
!!$         
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set the method-specific options.                              *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$
!!$         ! i) Set the maximum number of steps before restart (default=30).
!!$         
!!$         ! KSPGMRESSetRestart - Sets number of iterations at which GMRES,
!!$         !   FGMRES and LGMRES restarts.
!!$         !
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h"  
!!$         ! call KSPGMRESSetRestart(KSP ksp, PetscInt restart, &
!!$         !                         PetscErrorCode ierr)
!!$         !
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameters
!!$         !   ksp     - the Krylov space context
!!$         !   restart - integer restart value
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPGMRESSetRestart.html
!!$         ! or PETSc users manual, pp.65
!!$         
!!$         call KSPGMRESSetRestart(ksp, adjRestart, PETScIerr)
!!$         
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", "Error in KSPGMRESSetRestart")
!!$
!!$         ! ii) Set the iterative refinement type (default=none) for the
!!$         !     orthogonalization of the Hessenberg matrix
!!$         
!!$         ! KSPGMRESSetCGSRefinementType - Sets the type of iterative
!!$         !                refinement to use in the classical Gram Schmidt
!!$         !                orthogonalization. of the preconditioned problem.
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h"  
!!$         ! call KSPGMRESSetCGSRefinementType(KSP ksp,                     &
!!$         !                                KSPGMRESCGSRefinementType type, &
!!$         !                                PetscErrorCode ierr)
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameters
!!$         !   ksp  - the Krylov space context
!!$         !   type - the type of refinement: KSP_GMRES_CGS_REFINE_NEVER
!!$         !                                  KSP_GMRES_CGS_REFINE_IFNEEDED
!!$         !                                  KSP_GMRES_CGS_REFINE_ALWAYS
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPGMRESSetCGSRefinementType.html
!!$         ! or PETSc users manual, pp.65
!!$         
!!$         !      call KSPGMRESSetCGSRefinementType(ksp, &
!!$         !                     KSP_GMRES_CGS_REFINE_NEVER, PETScIerr)
!!$         
!!$         call KSPGMRESSetCGSRefinementType(ksp, &
!!$              KSP_GMRES_CGS_REFINE_IFNEEDED, PETScIerr)
!!$
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", &
!!$              "Error in KSPGMRESSetCGSRefinementType")
!!$
!!$
!!$
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set left- or right-side preconditioning (default=left)        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$
!!$         ! KSPSetPreconditionerSide - Sets the preconditioning side.
!!$         !
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPSetPreconditionerSide(KSP ksp,PCSide side,&
!!$         !                               PetscErrorCode ierr)
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameter
!!$         ! ksp  - iterative context obtained from KSPCreate()
!!$         !
!!$         ! Output Parameter
!!$         ! side - the preconditioning side, where side is one of
!!$         !          PC_LEFT - left preconditioning (default)
!!$         !          PC_RIGHT - right preconditioning
!!$         !          PC_SYMMETRIC - symmetric preconditioning
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPSetPreconditionerSide.html
!!$         ! or PETSc users manual, pp.66
!!$         select case(PCSide)
!!$         case(Right)
!!$            call KSPSetPreconditionerSide(ksp, PC_RIGHT, PETScIerr)
!!$         case(Left)
!!$            call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)
!!$         end select
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", &
!!$              "Error in KSPSetPreconditionerSide")
!!$
!!$      case(PETSCBICGStab)
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set the Krylov method.                                        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$      
!!$         ! KSPSetType - Builds KSP for a particular solver.
!!$         ! 
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPSetType(KSP ksp, KSPType type, PetscErrorCode ierr)
!!$         !
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameters
!!$         !   ksp  - the Krylov space context
!!$         !   type - a known method: KSPGMRES - GMRES (default)
!!$         !                          KSPCG    - Conjugate gradients
!!$         !                          KSPLSQR  -Least-squares method
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPSetType.html
!!$         ! or PETSc users manual, pp.65
!!$
!!$         call KSPSetType(ksp, KSPBCGS, PETScIerr)
!!$         
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", "Error in KSPSetType")
!!$   
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Get the Krylov method.                                        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$         ! KSPGetType - Gets the KSP type as a string from the KSP object.
!!$         ! 
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPGetType(KSP ksp, KSPType type, PetscErrorCode ierr)
!!$         !
!!$         ! Not collective
!!$         !
!!$         ! Input Parameters
!!$         !   ksp  - the Krylov space context
!!$         ! Output Parameters
!!$         !   type - KSP type
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPGetType.html
!!$         
!!$         if( PETScRank==0 ) then
!!$            call KSPGetType(ksp, kspType, PETScIerr)
!!$            if( PETScIerr/=0 ) &
!!$                 call terminate("createPETScKsp", "Error in KSPGetType")
!!$            write(*,10) kspType
!!$         endif
!!$         
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set left- or right-side preconditioning (default=left)        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$
!!$         ! KSPSetPreconditionerSide - Sets the preconditioning side.
!!$         !
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPSetPreconditionerSide(KSP ksp,PCSide side,&
!!$         !                               PetscErrorCode ierr)
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameter
!!$         ! ksp  - iterative context obtained from KSPCreate()
!!$         !
!!$         ! Output Parameter
!!$         ! side - the preconditioning side, where side is one of
!!$         !          PC_LEFT - left preconditioning (default)
!!$         !          PC_RIGHT - right preconditioning
!!$         !          PC_SYMMETRIC - symmetric preconditioning
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPSetPreconditionerSide.html
!!$         ! or PETSc users manual, pp.66
!!$         select case(PCSide)
!!$         case(Right)
!!$            call KSPSetPreconditionerSide(ksp, PC_RIGHT, PETScIerr)
!!$         case(Left)
!!$            call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)
!!$         end select
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", &
!!$              "Error in KSPSetPreconditionerSide")
!!$
!!$      case(PETSCCG)
!!$
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set the Krylov method.                                        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$      
!!$         ! KSPSetType - Builds KSP for a particular solver.
!!$         ! 
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPSetType(KSP ksp, KSPType type, PetscErrorCode ierr)
!!$         !
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameters
!!$         !   ksp  - the Krylov space context
!!$         !   type - a known method: KSPGMRES - GMRES (default)
!!$         !                          KSPCG    - Conjugate gradients
!!$         !                          KSPLSQR  -Least-squares method
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPSetType.html
!!$         ! or PETSc users manual, pp.65
!!$
!!$         call KSPSetType(ksp, KSPCG, PETScIerr)
!!$         
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", "Error in KSPSetType")
!!$   
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Get the Krylov method.                                        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$         ! KSPGetType - Gets the KSP type as a string from the KSP object.
!!$         ! 
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPGetType(KSP ksp, KSPType type, PetscErrorCode ierr)
!!$         !
!!$         ! Not collective
!!$         !
!!$         ! Input Parameters
!!$         !   ksp  - the Krylov space context
!!$         ! Output Parameters
!!$         !   type - KSP type
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPGetType.html
!!$         
!!$         if( PETScRank==0 ) then
!!$            call KSPGetType(ksp, kspType, PETScIerr)
!!$            if( PETScIerr/=0 ) &
!!$                 call terminate("createPETScKsp", "Error in KSPGetType")
!!$            write(*,10) kspType
!!$         endif
!!$         
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set left- or right-side preconditioning (default=left)        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$
!!$         ! KSPSetPreconditionerSide - Sets the preconditioning side.
!!$         !
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPSetPreconditionerSide(KSP ksp,PCSide side,&
!!$         !                               PetscErrorCode ierr)
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameter
!!$         ! ksp  - iterative context obtained from KSPCreate()
!!$         !
!!$         ! Output Parameter
!!$         ! side - the preconditioning side, where side is one of
!!$         !          PC_LEFT - left preconditioning (default)
!!$         !          PC_RIGHT - right preconditioning
!!$         !          PC_SYMMETRIC - symmetric preconditioning
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPSetPreconditionerSide.html
!!$         ! or PETSc users manual, pp.66
!!$         select case(PCSide)
!!$         case(Right)
!!$            call KSPSetPreconditionerSide(ksp, PC_RIGHT, PETScIerr)
!!$         case(Left)
!!$            call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)
!!$         end select
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", &
!!$              "Error in KSPSetPreconditionerSide")
!!$  
!!$         
!!$      case(PETSCFGMRES)
!!$
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set the Krylov method.                                        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$
!!$         ! KSPSetType - Builds KSP for a particular solver.
!!$         ! 
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPSetType(KSP ksp, KSPType type, PetscErrorCode ierr)
!!$         !
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameters
!!$         !   ksp  - the Krylov space context
!!$         !   type - a known method: KSPGMRES - GMRES (default)
!!$         !                          KSPCG    - Conjugate gradients
!!$         !                          KSPLSQR  -Least-squares method
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPSetType.html
!!$         ! or PETSc users manual, pp.65
!!$         
!!$         call KSPSetType(ksp, KSPFGMRES, PETScIerr)
!!$         
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", "Error in KSPSetType")  
!!$
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Get the Krylov method.                                        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$         ! KSPGetType - Gets the KSP type as a string from the KSP object.
!!$         ! 
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPGetType(KSP ksp, KSPType type, PetscErrorCode ierr)
!!$         !
!!$         ! Not collective
!!$         !
!!$         ! Input Parameters
!!$         !   ksp  - the Krylov space context
!!$         ! Output Parameters
!!$         !   type - KSP type
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPGetType.html
!!$         
!!$         if( PETScRank==0 ) then
!!$            call KSPGetType(ksp, kspType, PETScIerr)
!!$            if( PETScIerr/=0 ) &
!!$                 call terminate("createPETScKsp", "Error in KSPGetType")
!!$            write(*,10) kspType
!!$         endif
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set the method-specific options.                              *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$
!!$         ! i) Set the maximum number of steps before restart (default=30).
!!$         
!!$         ! KSPGMRESSetRestart - Sets number of iterations at which GMRES,
!!$         !   FGMRES and LGMRES restarts.
!!$         !
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h"  
!!$         ! call KSPGMRESSetRestart(KSP ksp, PetscInt restart, &
!!$         !                         PetscErrorCode ierr)
!!$         !
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameters
!!$         !   ksp     - the Krylov space context
!!$         !   restart - integer restart value
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPGMRESSetRestart.html
!!$         ! or PETSc users manual, pp.65
!!$         
!!$         call KSPGMRESSetRestart(ksp, adjRestart, PETScIerr)
!!$         
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", "Error in KSPGMRESSetRestart")
!!$
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set left- or right-side preconditioning (default=left)        *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$
!!$         ! KSPSetPreconditionerSide - Sets the preconditioning side.
!!$         !
!!$         ! Synopsis
!!$         !
!!$         ! #include "petscksp.h" 
!!$         ! call KSPSetPreconditionerSide(KSP ksp,PCSide side,&
!!$         !                               PetscErrorCode ierr)
!!$         ! Collective on KSP
!!$         !
!!$         ! Input Parameter
!!$         ! ksp  - iterative context obtained from KSPCreate()
!!$         !
!!$         ! Output Parameter
!!$         ! side - the preconditioning side, where side is one of
!!$         !          PC_LEFT - left preconditioning (default)
!!$         !          PC_RIGHT - right preconditioning
!!$         !          PC_SYMMETRIC - symmetric preconditioning
!!$         !
!!$         ! see .../petsc/docs/manualpages/KSP/KSPSetPreconditionerSide.html
!!$         ! or PETSc users manual, pp.66
!!$         select case(PCSide)
!!$         case(Right)
!!$            call KSPSetPreconditionerSide(ksp, PC_RIGHT, PETScIerr)
!!$         case(Left)
!!$            call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)
!!$         end select
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", &
!!$              "Error in KSPSetPreconditionerSide")
!!$
!!$      end select
!!$
!!$!
!!$!     *****************************************************************
!!$!     *                                                               *
!!$!     * Set some options that are not solver specific                 *
!!$!     *                                                               *
!!$!     *****************************************************************
!!$!
!!$!
!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     * Set the convergence tolerances.                                *
!!$!     *                                                                *
!!$!     * rtol = decrease of the residual norm relative to the norm of   *
!!$!     *        the right-hand side                                     *
!!$!     *        (setting to PETSC_DEFAULT_DOUBLE_PRECISION = 10e-5)     *
!!$!     * atol = absolute size of the residual norm                      *
!!$!     *        (setting to PETSC_DEFAULT_DOUBLE_PRECISION = 10e-50)    *
!!$!     *                                                                *
!!$!     * dtol = relative increase in the residual                       *
!!$!     *        (setting to PETSC_DEFAULT_DOUBLE_PRECISION = 10e+5)     *
!!$!     *                                                                *
!!$!     * maxIter = maximum number of allowable iterations               *
!!$!     *           (setting to PETSC_DEFAULT_INTEGER = 10e+4)           *
!!$!     *                                                                *
!!$!     ******************************************************************
!!$!
!!$      ! KSPSetTolerances - Sets the relative, absolute, divergence, and
!!$      !                    maximum iteration tolerances used by the
!!$      !                    default KSP convergence testers.
!!$      ! Synopsis
!!$      !
!!$      ! #include "petscksp.h" 
!!$      ! call KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,&
!!$      !              PetscReal dtol,PetscInt maxits,PetscErrorCode ierr)
!!$      !
!!$      ! Collective on KSP
!!$      !
!!$      ! Input Parameters
!!$      !   ksp    - the Krylov subspace context
!!$      !   rtol   - the relative convergence tolerance (relative decrease
!!$      !            in the residual norm)
!!$      !   abstol - the absolute convergence tolerance (absolute size of
!!$      !            the residual norm)
!!$      !   dtol   - the divergence tolerance (amount residual can
!!$      !            increase before KSPDefaultConverged() concludes that
!!$      !            the method is diverging)
!!$      !   maxits - maximum number of iterations to use
!!$      !
!!$      ! see .../petsc/docs/manualpages/KSP/KSPSetTolerances.html
!!$      ! or PETSc users manual, pp.67
!!$
!!$      call KSPSetTolerances(ksp, adjRelTol, adjAbsTol, adjDivTol, &
!!$                            adjMaxIter, PETScIerr)
!!$
!!$      if( PETScIerr/=0 ) &
!!$        call terminate("createPETScKsp", "Error in KSPSetTolerances")
!!$
!!$      ! KSPGetTolerances - Gets the relative, absolute, divergence, and
!!$      !                    maximum iteration tolerances used by the
!!$      !                    default KSP convergence tests.
!!$      ! Synopsis
!!$      !
!!$      ! #include "petscksp.h" 
!!$      ! call KSPGetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,&
!!$      !           PetscReal *dtol,PetscInt *maxits, PetscErrorCode ierr)
!!$      !
!!$      ! Not Collective
!!$      !
!!$      ! Input Parameter
!!$      !   ksp -the Krylov subspace context
!!$      ! Output Parameters
!!$      !   rtol   - the relative convergence tolerance
!!$      !   abstol - the absolute convergence tolerance
!!$      !   dtol 	 - the divergence tolerance
!!$      !   maxits - maximum number of iterations
!!$      !
!!$      ! Notes
!!$      ! The user can specify PETSC_NULL for
!!$      !   any parameter that is not needed.
!!$      !
!!$      ! see .../petsc/docs/manualpages/KSP/KSPGetTolerances.html
!!$
!!$      if( PETScRank==0 ) then
!!$         call KSPGetTolerances(ksp, rTol, aTol, dTol, mIts, PETScIerr)
!!$         if( PETScIerr/=0 ) &
!!$              call terminate("createPETScKsp", "Error in KSPGetTolerances")
!!$         write(*,20) rTol, aTol, dTol
!!$         write(*,21) mIts
!!$      endif
!!$
!!$!
!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     * Set the preconditioner. The PETSc libraries supports the       *
!!$!     * following preconditioners:                                     *
!!$!     *      PCJACOBI    - Jabobi                                      *
!!$!     *      PCBJACOBI   - Block Jacobi                                *
!!$!     *      PCSOR       - SOR                                         *
!!$!     *      PCEISENSTAT - SOR with Eisenstat trick                    *
!!$!     *      PCICC       - incomplete Cholesky                         *
!!$!     *      PCILU       - incomplete LU                               *
!!$!     *      PCASM       - additive Schwartz                           *
!!$!     *      PCKSP       - linear solver                               *
!!$!     *      PCCOMPOSITE - combination of preconditioners              *
!!$!     *      PCLU        - LU                                          *
!!$!     *      PCCHOLESKY  - Cholesky                                    *
!!$!     *      PCNONE      - no preconditioning                          *
!!$!     *      PCMG        - multigrid preconditioning                   *
!!$!     *                                                                *
!!$!     * For more information, see PETSc users manual, pp.70-77.        *
!!$!     *                                                                *
!!$!     ****************************************************************** 
!!$
!!$      ! By extracting the KSP and PC contexts from the KSP context,
!!$      ! we can then directly directly call any KSP and PC routines
!!$      ! to set various options (PETSc users manual, pp.64).
!!$      
!!$      ! KSPGetPC - Returns a pointer to the preconditioner context set
!!$      !            with KSPSetPC().
!!$      ! Synopsis
!!$      !
!!$      ! #include "petscksp.h" 
!!$      ! call KSPGetPC(KSP ksp,PC *pc,PetscErrorCode ierr)
!!$      !
!!$      ! Not Collective
!!$      !
!!$      ! Input Parameters
!!$      !   ksp - iterative context obtained from KSPCreate()
!!$      !
!!$      ! Output Parameter
!!$      !   pc - preconditioner context
!!$      !
!!$      ! see .../petsc/docs/manualpages/KSP/KSPGetPC.html
!!$      
!!$      call KSPGetPC(ksp, pc, PETScIerr)
!!$      
!!$      if( PETScIerr/=0 ) &
!!$           call terminate("createPETScKsp", "Error in KSPGetPC")
!!$
!!$!
!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     * Select the preconditioning method.                             *
!!$!     *                                                                *
!!$!     ****************************************************************** 
!!$!
!!$
!!$select case(PreCondType)
!!$   
!!$case(Jacobi)
!!$
!!$   !call terminate("createPETScKsp", "Jacobi Preconditioner Not implemented")
!!$   ! PCSetType - Builds PC for a particular preconditioner.
!!$   !
!!$   ! Synopsis
!!$   !
!!$   ! #include "petscpc.h" 
!!$   ! call PCSetType(PC pc, PCType type, PetscErrorCode ierr)
!!$   !
!!$   ! Collective on PC
!!$   !
!!$   ! Input Parameter
!!$   !   pc   - the preconditioner context.
!!$   !   type - a known method 
!!$   !
!!$   ! see .../petsc/docs/manualpages/PC/PCSetType.html
!!$   
!!$   call PCSetType( pc, PCJACOBI, PETScIerr)
!!$   
!!$   ! set tolerances on subksp (is this really needed???)
!!$   call KSPSetTolerances(ksp, 1.e-8, adjAbsTol, adjDivTol, &
!!$        adjMaxIter, PETScIerr)
!!$   
!!$
!!$   
!!$   !set Scaling Factor Type
!!$   select case(ScaleType)
!!$   case(Normal)
!!$      continue
!!$   case(RowMax)
!!$      call PCJacobiSetUseRowMax(pc, PETScIerr)
!!$   case(RowSum)
!!$      call PCJacobiSetUseRowSum(pc, PETScIerr)
!!$   case(RowAbs)
!!$      call PCJacobiSetUseAbs(pc, PETScIerr)
!!$   end select

!!$   !set matrix ordering
!!$   
!!$   select case(MatrixOrdering)
!!$   case(Natural)
!!$      call PCFactorSetMatOrderingtype( pc, MATORDERING_NATURAL, PETScIerr )
!!$   case(ReverseCuthillMckee)
!!$      call PCFactorSetMatOrderingtype( pc, MATORDERING_RCM, PETScIerr )
!!$   case(NestedDissection)
!!$      call PCFactorSetMatOrderingtype( pc, MATORDERING_ND, PETScIerr )
!!$   case(OnewayDissection )
!!$      call PCFactorSetMatOrderingtype( pc, MATORDERING_1WD, PETScIerr )
!!$   case( QuotientMinimumDegree)
!!$      call PCFactorSetMatOrderingtype( pc, MATORDERING_QMD, PETScIerr )
!!$   end select
!!$
!!$   !Set the iteration Monitor
!!$   if (setMonitor) then
!!$      call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
!!$           PETSC_NULL_FUNCTION, PETScIerr)
!!$      !call KSPMonitorSet(ksp,MyKSPMonitor, PETScIerr)
!!$   endif
!!$
!!$case(BlockJacobi)
!!$   ! PCSetType - Builds PC for a particular preconditioner.
!!$   !
!!$   ! Synopsis
!!$   !
!!$   ! #include "petscpc.h" 
!!$   ! call PCSetType(PC pc, PCType type, PetscErrorCode ierr)
!!$   !
!!$   ! Collective on PC
!!$   !
!!$   ! Input Parameter
!!$   !   pc   - the preconditioner context.
!!$   !   type - a known method 
!!$   !
!!$   ! see .../petsc/docs/manualpages/PC/PCSetType.html
!!$
!!$   call PCSetType( pc, PCBJACOBI, PETScIerr)
!!$   call PCSetUp(pc,PETScIerr)
!!$   Nsub = 1!_intType
!!$   length = nCellsLocal*nw
!!$   call PCBJacobiSetLocalBlocks(pc,Nsub,length,PETScIerr)
!!$
!!$   !Setup KSP context before calling local subdomain ksp contexts
!!$   call KSPSetUp(ksp, PETScIerr)
!!$
!!$   !Setup local subcontexts
!!$   call PCBJacobiGetSubKSP(pc,nlocal,first,subksp,PETScIerr)
!!$   
!!$   !Setup Local ILU precondtioner
!!$   call KSPGetPC( subksp, subpc, PETScIerr )
!!$   
!!$   select case(LocalPCType)
!!$   case(ILU)
!!$      call PCSetType( subpc, PCILU, PETScIerr )
!!$   case(ICC)
!!$      call PCSetType( subpc, PCICC, PETScIerr )
!!$   case(LU)
!!$      call PCSetType( subpc, PCLU, PETScIerr )
!!$   case(Cholesky)
!!$      call PCSetType( subpc, PCCHOLESKY, PETScIerr )
!!$   end select
!!$
!!$!set matrix ordering
!!$
!!$   select case(MatrixOrdering)
!!$   case(Natural)
!!$      call PCFactorSetMatOrderingtype( subpc, MATORDERING_NATURAL, PETScIerr )
!!$   case(ReverseCuthillMckee)
!!$      call PCFactorSetMatOrderingtype( subpc, MATORDERING_RCM, PETScIerr )
!!$   case(NestedDissection)
!!$      call PCFactorSetMatOrderingtype( subpc, MATORDERING_ND, PETScIerr )
!!$   case(OnewayDissection )
!!$      call PCFactorSetMatOrderingtype( subpc, MATORDERING_1WD, PETScIerr )
!!$   case( QuotientMinimumDegree)
!!$      call PCFactorSetMatOrderingtype( subpc, MATORDERING_QMD, PETScIerr )
!!$   end select
!!$
!!$   !Set ILU parameters
!!$   call PCFactorSetLevels( subpc, fillLevel , PETScIerr)!  set 1 level of fill
!!$   !call PCFactorSetFill( subpc, 3, PETScIerr )
!!$
!!$   !Set local contexts to preconditioner's only
!!$   call KSPSetType(subksp, KSPPREONLY, PETScIerr)
!!$   
!!$   ! set tolerances on subksp (is this really needed???)
!!$   call KSPSetTolerances(subksp, 1.e-8, adjAbsTol, adjDivTol, &
!!$        adjMaxIter, PETScIerr)
!!$   
!!$   if(setMonitor)then
!!$      !Set the convergence monitors
!!$      call KSPMonitorSet(subksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
!!$           PETSC_NULL_FUNCTION, PETScIerr)
!!$      
!!$      call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
!!$           PETSC_NULL_FUNCTION, PETScIerr)
!!$   endif
!!$
!!$case(AdditiveSchwartz)
!!$
!!$   !call terminate("createPETScKsp", "ASM Preconditioner Not implemented")
!!$      ! PCSetType - Builds PC for a particular preconditioner.
!!$      !
!!$      ! Synopsis
!!$      !
!!$      ! #include "petscpc.h" 
!!$      ! call PCSetType(PC pc, PCType type, PetscErrorCode ierr)
!!$      !
!!$      ! Collective on PC
!!$      !
!!$      ! Input Parameter
!!$      !   pc   - the preconditioner context.
!!$      !   type - a known method 
!!$      !
!!$      ! see .../petsc/docs/manualpages/PC/PCSetType.html
!!$      call PCSetType( pc, PCASM, PETScIerr)
!!$
!!$      !setup a basic overlaping scheme, more detailed scheme to be used later
!!$      
!!$      call PCASMSetOverlap(pc,overlap,PETScIerr)
!!$      
!!$      !Setup KSP context before calling local subdomain ksp contexts
!!$      call KSPSetUp(ksp, PETScIerr);
!!$      
!!$      !get the subdomain contexts
!!$      call PCASMGetSubKSP( pc, nlocal,  first, subksp, PETScIerr )
!!$      
!!$      !Setup the subdomain preconditioner
!!$      call KSPGetPC( subksp, subpc, PETScIerr )
!!$      !Set subdomain preconditioner type
!!$      
!!$      select case(LocalPCType)
!!$      case(ILU)
!!$         call PCSetType( subpc, PCILU, PETScIerr )
!!$      case(ICC)
!!$         call PCSetType( subpc, PCICC, PETScIerr )
!!$      case(LU)
!!$         call PCSetType( subpc, PCLU, PETScIerr )
!!$      case(Cholesky)
!!$         call PCSetType( subpc, PCCHOLESKY, PETScIerr )
!!$      end select
!!$      
!!$      !set matrix ordering
!!$      
!!$      select case(MatrixOrdering)
!!$      case(Natural)
!!$         call PCFactorSetMatOrderingtype( subpc, MATORDERING_NATURAL, PETScIerr )
!!$      case(ReverseCuthillMckee)
!!$         call PCFactorSetMatOrderingtype( subpc, MATORDERING_RCM, PETScIerr )
!!$      case(NestedDissection)
!!$         call PCFactorSetMatOrderingtype( subpc, MATORDERING_ND, PETScIerr )
!!$      case(OnewayDissection )
!!$         call PCFactorSetMatOrderingtype( subpc, MATORDERING_1WD, PETScIerr )
!!$      case( QuotientMinimumDegree)
!!$         call PCFactorSetMatOrderingtype( subpc, MATORDERING_QMD, PETScIerr )
!!$      end select
!!$
!!$      !Set ILU parameters
!!$      call PCFactorSetLevels( subpc, fillLevel , PETScIerr)!  set 1 level of fill
!!$      !call PCFactorSetFill( subpc, 3, PETScIerr )
!!$      
!!$      !Set local contexts to preconditioner's only
!!$      call KSPSetType(subksp, KSPPREONLY, PETScIerr)
!!$
!!$      ! set tolerances on subksp (is this really needed???)
!!$      call KSPSetTolerances(subksp, 1.e-8, adjAbsTol, adjDivTol, &
!!$           adjMaxIter, PETScIerr)
!!$      
!!$      if(setMonitor)then
!!$         !Set the convergence monitors
!!$!         call KSPMonitorSet(subksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
!!$!              PETSC_NULL_FUNCTION, PETScIerr)
!!$         
!!$         call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
!!$              PETSC_NULL_FUNCTION, PETScIerr)
!!$
!!$!         call KSPMonitorSet(ksp,KSPMonitorSingularValue,PETSC_NULL_OBJECT,&
!!$!              PETSC_NULL_FUNCTION, PETScIerr)
!!$
!!$      endif
!!$end select
!!$
!!$!
!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     * Get the preconditioning method details and print to screen.    *
!!$!     *                                                                *
!!$!     ****************************************************************** 
!!$!
!!$!
!!$
!!$! PCGetType - Gets the PC method type and name (as a string)
!!$!             from the PC context.
!!$! Synopsis
!!$!
!!$! #include "petscpc.h" 
!!$! call PCGetType(PC pc,PCType meth, PetscErrorCode ierr)
!!$!
!!$! Not Collective
!!$!
!!$! Input Parameter
!!$!   pc   - the preconditioner context
!!$!
!!$! Output Parameter
!!$!   meth - name of preconditioner 
!!$!
!!$! see .../petsc/docs/manualpages/PC/PCGetType.html
!!$
!!$if( PETScRank==0 ) then
!!$   call PCGetType(pc, pcType, PETScIerr)
!!$   if( PETScIerr/=0 ) &
!!$        call terminate("createPETScKsp", "Error in PCGetType")
!!$   write(*,30) pcType
!!$endif

! Output formats.

10 format("# ... KSP properties:",/, &
        "#",7x,"type        :",1x,a)
20 format("#",7x,"tolerances  :",1x,"rel =",1x,es7.1,/, &
        "#",21x,"abs =",1x,es7.1,/, &
        "#",21x,"div =",1x,es7.1)
21 format("#",7x,"max.iter.   :",1x,i4)
30 format("#",7x,"precond.type:",1x,a)
40 format(a,1x,i3,a,1x,i6,a,1x,i6,1x,a,1x,i6) ! ownership
#endif

end subroutine createPETScKsp
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine MyKSPMonitor(myKsp, n, rnorm, dummy, ierr)
!!$  !
!!$  !     ******************************************************************
!!$  !     *                                                                *
!!$  !     * This is a user-defined routine for monitoring the KSP          *
!!$  !     * iterative solvers. Instead of outputing the L2-norm at every   *
!!$  !     * iteration (default PETSc monitor), it only does it every       *
!!$  !     * 'adjMonStep' iterations.                                       *
!!$  !     *                                                                *
!!$  !     ******************************************************************
!!$  !
!!$  use ADjointPETSc
!!$  use inputADjoint
!!$  implicit none
!!$  !
!!$  !     Subroutine arguments.
!!$  !
!!$  ! myKsp - Iterative context
!!$  ! n     - Iteration number
!!$  ! rnorm - 2-norm (preconditioned) residual value
!!$  ! dummy - Optional user-defined monitor context (unused here)
!!$  ! ierr  - Return error code
!!$  
!!$  real(kind=realType), pointer, dimension(:,:) :: myKsp
!!$  integer(kind=intType) :: n, dummy, ierr
!!$  real(kind=realType)   :: rnorm
!!$  !
!!$  !     ******************************************************************
!!$  !     *                                                                *
!!$  !     * Begin execution.                                               *
!!$  !     *                                                                *
!!$  !     ******************************************************************
!!$  !
!!$#ifndef USE_NO_PETSC
!!$  
!!$  ! Write the residual norm to stdout every adjMonStep iterations.
!!$  
!!$  if( mod(n,adjMonStep)==0 ) then
!!$        if( PETScRank==0 ) write(*,10) n, rnorm
!!$     end if
!!$
!!$     ierr = 0
!!$     
!!$     ! Output format.
!!$     
!!$10   format(i4,1x,'KSP Residual norm',1x,e10.4)
!!$     
!!$#endif
!!$     
!!$   end subroutine MyKSPMonitor
   
