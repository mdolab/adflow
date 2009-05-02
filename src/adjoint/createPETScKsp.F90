!
!     ******************************************************************
!     *                                                                *
!     * File:          createPETScKsp.F90                              *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 12-15-2005                                      *
!     * Last modified: 02-07-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine createPETScKsp
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
!      use internalCommunication
      use blockPointers
!      use indices         ! nw
      use flowVarRefState  !nw
      implicit none

!#include "include/finclude/petsc.h"
!#include "include/finclude/petsc.h"
!#include "include/finclude/petscis.h"

!
!     User-defined Fortran routine.
!
      external MyKSPMonitor
!
!     Local variables.
!
      integer(kind=intType) :: level = 1	
      character(len=10)     :: kspType, pcType
      real(kind=realType)   :: rTol, aTol, dTol
      integer(kind=intType) :: mIts
!      integer(kind=intType) :: i!,nlocal,
      integer       :: nn, iLow, iHigh!,length,overlap=1
      integer(kind=intType),dimension(:),allocatable::localASMRange
      integer(kind=intType) ::nPenaltyLocal

      integer(kind=intType) :: ii, sps=1, nNode,counter,i,j,k
      integer(kind=intType) :: b1, b2, s1, s2, i1, i2, j1, j2, ir
      integer(kind=intType), dimension(:),   pointer :: block1, block2
      integer(kind=intType), dimension(:),   pointer :: sub1,   sub2
      integer(kind=intType), dimension(:),   pointer :: indPen
      integer(kind=intType), dimension(:,:), pointer :: ind1,   ind2

      logical :: ASM =.False.!.True.!.False.
      logical :: BJACOBI =.True.!.False.!

      logical :: GMRES =.True.!.False.
      logical :: FGMRES =.False.!.True.!.False.
      logical :: BICGSTAB =.False.!.True.!.False.

!      PetscInt  ierr 
!      PetscInt  indices2(5),rank,n
!      integer, pointer :: idx(:)

!       PetscErrorCode ierr
!       PetscInt indices2(5),n,index1,index5
!       PetscMPIInt rank
!       PetscOffset ix
!       IS          is	
      	
       integer ierr
       integer indices2(5),n,index1,index5
       integer  rank
!       PetscOffset ix
!       IS          is	

!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! PETSc macros are lost and have to be redefined.
      ! They were extracted from: <PETSc root dir>/include/petscksp.h
      !                                                   /petscpc.h
#define KSPGMRES      "gmres"   
#define PCILU         "ilu"
#define MATORDERING_RCM       "rcm"
#define PCASM        "asm"
#define PCBJACOBI       "bjacobi"
#define KSPCG         "cg"
#define PETSC_NULL           0
#define PETSC_DEFAULT        -2
#define KSPPREONLY    "preonly"
#define KSPBCGS       "bcgs"
#define PCLU 'lu'
#define KSPFGMRES     "fgmres"
#define MATORDERING_NATURAL       "natural"
#define MATORDERING_ND        "nd"
#define MATORDERING_1WD       "1wd"
#define MATORDERING_RCM       "rcm"
#define MATORDERING_QMD       "qmd"

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

!
!     ******************************************************************
!     *                                                                *
!     * Solve the linear system of equations dRdWT . psi = dIdW using  *
!     * preconditioned GMRES.                                          *
!     *                                                                *
!     ******************************************************************
!
      ! KSPSetOperators - Sets the matrix associated with the linear
      !                   system and a (possibly) different one
      !                   associated with the preconditioner.
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetOperators(KSP ksp,Mat Amat,Mat Pmat, &
      !                      MatStructure flag,PetscErrorCode ierr)
      !
      ! Collective on KSP and Mat
      !
      ! Input Parameters
      !   ksp  - the KSP context
      !   Amat - the matrix associated with the linear system
      !   Pmat - the matrix to be used in constructing the 
      !          preconditioner, usually the same as Amat.
      !   flag - flag indicating information about the preconditioner
      !          matrix structure during successive linear solves.
      !          This flag is ignored the first time a linear system is
      !          solved, and thus is irrelevant when solving just one
      !          linear system.
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetOperators.html
      ! or PETSc users manual, pp.63

      ! Here the matrix that defines the linear system
      ! also serves as the preconditioning matrix.

      call KSPSetOperators(ksp,dRdW,dRdW, &
                           DIFFERENT_NONZERO_PATTERN,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("solveADjointPETSc", "Error in KSPSetOperators.")

!
!     ******************************************************************
!     *                                                                *
!     * Set the various options for the KSP.                           *
!     *                                                                *
!     ******************************************************************
!
      ! KSPSetFromOptions - Sets KSP options from the options database.
      !   This routine must be called before KSPSetUp() if the user is
      !   to be allowed to set the Krylov type.
      !
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetFromOptions(KSP ksp, PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp - the Krylov space context
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetFromOptions.html

      call KSPSetFromOptions(ksp, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPSetFromOptions")

      ! >>> The following statements are OPTIONAL

!***********************************
! GMRES   
!************************************
      if(GMRES)then
       ! Set the Krylov method.

      ! KSPSetType - Builds KSP for a particular solver.
      ! 
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetType(KSP ksp, KSPType type, PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp  - the Krylov space context
      !   type - a known method: KSPGMRES - GMRES (default)
      !                          KSPCG    - Conjugate gradients
      !                          KSPLSQR  -Least-squares method
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetType.html
      ! or PETSc users manual, pp.65

      call KSPSetType(ksp, KSPGMRES, PETScIerr)
      

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPSetType")

      ! Get the Krylov method.

      ! KSPGetType - Gets the KSP type as a string from the KSP object.
      ! 
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPGetType(KSP ksp, KSPType type, PetscErrorCode ierr)
      !
      ! Not collective
      !
      ! Input Parameters
      !   ksp  - the Krylov space context
      ! Output Parameters
      !   type - KSP type
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGetType.html

      if( PETScRank==0 ) then
        call KSPGetType(ksp, kspType, PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("createPETScKsp", "Error in KSPGetType")
        write(*,10) kspType
      endif

      ! Set the method-specific options.

      ! i) Set the maximum number of steps before restart (default=30).

      ! KSPGMRESSetRestart - Sets number of iterations at which GMRES,
      !   FGMRES and LGMRES restarts.
      !
      ! Synopsis
      !
      ! #include "petscksp.h"  
      ! call KSPGMRESSetRestart(KSP ksp, PetscInt restart, &
      !                         PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp     - the Krylov space context
      !   restart - integer restart value
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGMRESSetRestart.html
      ! or PETSc users manual, pp.65

      call KSPGMRESSetRestart(ksp, adjRestart, PETScIerr)
      !call KSPGMRESSetRestart(ksp, 500, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPGMRESSetRestart")



      ! ii) Set the iterative refinement type (default=none) for the
      !     orthogonalization of the Hessenberg matrix

      ! KSPGMRESSetCGSRefinementType - Sets the type of iterative
      !                refinement to use in the classical Gram Schmidt
      !                orthogonalization. of the preconditioned problem.
      ! Synopsis
      !
      ! #include "petscksp.h"  
      ! call KSPGMRESSetCGSRefinementType(KSP ksp,                     &
      !                                KSPGMRESCGSRefinementType type, &
      !                                PetscErrorCode ierr)
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp  - the Krylov space context
      !   type - the type of refinement: KSP_GMRES_CGS_REFINE_NEVER
      !                                  KSP_GMRES_CGS_REFINE_IFNEEDED
      !                                  KSP_GMRES_CGS_REFINE_ALWAYS
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGMRESSetCGSRefinementType.html
      ! or PETSc users manual, pp.65

!      call KSPGMRESSetCGSRefinementType(ksp, &
!                     KSP_GMRES_CGS_REFINE_NEVER, PETScIerr)

      call KSPGMRESSetCGSRefinementType(ksp, &
                     KSP_GMRES_CGS_REFINE_IFNEEDED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", &
                       "Error in KSPGMRESSetCGSRefinementType")

      ! Set left- or right-side preconditioning (default=left)

      ! KSPSetPreconditionerSide - Sets the preconditioning side.
      !
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetPreconditionerSide(KSP ksp,PCSide side,&
      !                               PetscErrorCode ierr)
      ! Collective on KSP
      !
      ! Input Parameter
      ! ksp  - iterative context obtained from KSPCreate()
      !
      ! Output Parameter
      ! side - the preconditioning side, where side is one of
      !          PC_LEFT - left preconditioning (default)
      !          PC_RIGHT - right preconditioning
      !          PC_SYMMETRIC - symmetric preconditioning
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetPreconditionerSide.html
      ! or PETSc users manual, pp.66

      call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", &
                       "Error in KSPSetPreconditionerSide")


      elseif(FGMRES)then

!*************************
!FGMRES
!*************************
       ! Set the Krylov method.

      ! KSPSetType - Builds KSP for a particular solver.
      ! 
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetType(KSP ksp, KSPType type, PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp  - the Krylov space context
      !   type - a known method: KSPGMRES - GMRES (default)
      !                          KSPCG    - Conjugate gradients
      !                          KSPLSQR  -Least-squares method
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetType.html
      ! or PETSc users manual, pp.65

      call KSPSetType(ksp, KSPFGMRES, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPSetType")

      ! Get the Krylov method.

      ! KSPGetType - Gets the KSP type as a string from the KSP object.
      ! 
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPGetType(KSP ksp, KSPType type, PetscErrorCode ierr)
      !
      ! Not collective
      !
      ! Input Parameters
      !   ksp  - the Krylov space context
      ! Output Parameters
      !   type - KSP type
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGetType.html

      if( PETScRank==0 ) then
        call KSPGetType(ksp, kspType, PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("createPETScKsp", "Error in KSPGetType")
        write(*,10) kspType
      endif

      ! Set the method-specific options.

      ! i) Set the maximum number of steps before restart (default=30).

      ! KSPGMRESSetRestart - Sets number of iterations at which GMRES,
      !   FGMRES and LGMRES restarts.
      !
      ! Synopsis
      !
      ! #include "petscksp.h"  
      ! call KSPGMRESSetRestart(KSP ksp, PetscInt restart, &
      !                         PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp     - the Krylov space context
      !   restart - integer restart value
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGMRESSetRestart.html
      ! or PETSc users manual, pp.65

      !call KSPGMRESSetRestart(ksp, adjRestart, PETScIerr)
      call KSPGMRESSetRestart(ksp, 500, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPGMRESSetRestart")



      ! ii) Set the iterative refinement type (default=none) for the
      !     orthogonalization of the Hessenberg matrix

      ! KSPGMRESSetCGSRefinementType - Sets the type of iterative
      !                refinement to use in the classical Gram Schmidt
      !                orthogonalization. of the preconditioned problem.
      ! Synopsis
      !
      ! #include "petscksp.h"  
      ! call KSPGMRESSetCGSRefinementType(KSP ksp,                     &
      !                                KSPGMRESCGSRefinementType type, &
      !                                PetscErrorCode ierr)
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp  - the Krylov space context
      !   type - the type of refinement: KSP_GMRES_CGS_REFINE_NEVER
      !                                  KSP_GMRES_CGS_REFINE_IFNEEDED
      !                                  KSP_GMRES_CGS_REFINE_ALWAYS
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGMRESSetCGSRefinementType.html
      ! or PETSc users manual, pp.65

!      call KSPGMRESSetCGSRefinementType(ksp, &
!                     KSP_GMRES_CGS_REFINE_NEVER, PETScIerr)

!      call KSPGMRESSetCGSRefinementType(ksp, &
!                     KSP_GMRES_CGS_REFINE_IFNEEDED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", &
                       "Error in KSPGMRESSetCGSRefinementType")

      ! Set left- or right-side preconditioning (default=left)

      ! KSPSetPreconditionerSide - Sets the preconditioning side.
      !
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetPreconditionerSide(KSP ksp,PCSide side,&
      !                               PetscErrorCode ierr)
      ! Collective on KSP
      !
      ! Input Parameter
      ! ksp  - iterative context obtained from KSPCreate()
      !
      ! Output Parameter
      ! side - the preconditioning side, where side is one of
      !          PC_LEFT - left preconditioning (default)
      !          PC_RIGHT - right preconditioning
      !          PC_SYMMETRIC - symmetric preconditioning
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetPreconditionerSide.html
      ! or PETSc users manual, pp.66

      call KSPSetPreconditionerSide(ksp, PC_RIGHT, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", &
                       "Error in KSPSetPreconditionerSide")

     
     
      elseif(BICGSTAB) then
!*******************************
!BiCGStab
!*******************************

      ! Set the Krylov method.

      ! KSPSetType - Builds KSP for a particular solver.
      ! 
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetType(KSP ksp, KSPType type, PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp  - the Krylov space context
      !   type - a known method: KSPGMRES - GMRES (default)
      !                          KSPCG    - Conjugate gradients
      !                          KSPLSQR  -Least-squares method
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetType.html
      ! or PETSc users manual, pp.65

      call KSPSetType(ksp, KSPBCGS, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPSetType")

      ! Get the Krylov method.

      ! KSPGetType - Gets the KSP type as a string from the KSP object.
      ! 
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPGetType(KSP ksp, KSPType type, PetscErrorCode ierr)
      !
      ! Not collective
      !
      ! Input Parameters
      !   ksp  - the Krylov space context
      ! Output Parameters
      !   type - KSP type
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGetType.html

      if( PETScRank==0 ) then
        call KSPGetType(ksp, kspType, PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("createPETScKsp", "Error in KSPGetType")
        write(*,10) kspType
      endif


      ! Set left- or right-side preconditioning (default=left)

      ! KSPSetPreconditionerSide - Sets the preconditioning side.
      !
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetPreconditionerSide(KSP ksp,PCSide side,&
      !                               PetscErrorCode ierr)
      ! Collective on KSP
      !
      ! Input Parameter
      ! ksp  - iterative context obtained from KSPCreate()
      !
      ! Output Parameter
      ! side - the preconditioning side, where side is one of
      !          PC_LEFT - left preconditioning (default)
      !          PC_RIGHT - right preconditioning
      !          PC_SYMMETRIC - symmetric preconditioning
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetPreconditionerSide.html
      ! or PETSc users manual, pp.66

      call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", &
                       "Error in KSPSetPreconditionerSide")



      else

      call terminate("createPETScKsp", &
                       "No valid solver specified")
 
      endif
!*********************************
!end solver select
!*********************************

!
!     ******************************************************************
!     *                                                                *
!     * Set the convergence tolerances.                                *
!     *                                                                *
!     * rtol = decrease of the residual norm relative to the norm of   *
!     *        the right-hand side                                     *
!     *        (setting to PETSC_DEFAULT_DOUBLE_PRECISION = 10e-5)     *
!     * atol = absolute size of the residual norm                      *
!     *        (setting to PETSC_DEFAULT_DOUBLE_PRECISION = 10e-50)    *
!     *                                                                *
!     * dtol = relative increase in the residual                       *
!     *        (setting to PETSC_DEFAULT_DOUBLE_PRECISION = 10e+5)     *
!     *                                                                *
!     * maxIter = maximum number of allowable iterations               *
!     *           (setting to PETSC_DEFAULT_INTEGER = 10e+4)           *
!     *                                                                *
!     ******************************************************************
!
      ! KSPSetTolerances - Sets the relative, absolute, divergence, and
      !                    maximum iteration tolerances used by the
      !                    default KSP convergence testers.
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,&
      !              PetscReal dtol,PetscInt maxits,PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp    - the Krylov subspace context
      !   rtol   - the relative convergence tolerance (relative decrease
      !            in the residual norm)
      !   abstol - the absolute convergence tolerance (absolute size of
      !            the residual norm)
      !   dtol   - the divergence tolerance (amount residual can
      !            increase before KSPDefaultConverged() concludes that
      !            the method is diverging)
      !   maxits - maximum number of iterations to use
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetTolerances.html
      ! or PETSc users manual, pp.67

      call KSPSetTolerances(ksp, adjRelTol, adjAbsTol, adjDivTol, &
                            2*adjMaxIter, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPSetTolerances")

      ! KSPGetTolerances - Gets the relative, absolute, divergence, and
      !                    maximum iteration tolerances used by the
      !                    default KSP convergence tests.
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPGetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,&
      !           PetscReal *dtol,PetscInt *maxits, PetscErrorCode ierr)
      !
      ! Not Collective
      !
      ! Input Parameter
      !   ksp -the Krylov subspace context
      ! Output Parameters
      !   rtol   - the relative convergence tolerance
      !   abstol - the absolute convergence tolerance
      !   dtol 	 - the divergence tolerance
      !   maxits - maximum number of iterations
      !
      ! Notes
      ! The user can specify PETSC_NULL for
      !   any parameter that is not needed.
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGetTolerances.html

      if( PETScRank==0 ) then
        call KSPGetTolerances(ksp, rTol, aTol, dTol, mIts, PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("createPETScKsp", "Error in KSPGetTolerances")
        write(*,20) rTol, aTol, dTol
        write(*,21) mIts
      endif


      ! Convergence monitoring (default = silent mode).

      ! KSPSetMonitor - Sets an ADDITIONAL function to be called at
      !               every iteration to monitor the residual/error etc.
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetMonitor(KSP ksp, &
      !        PetscErrorCode (*monitor)(KSP,PetscInt,PetscReal,void*),&
      !        void *mctx,PetscErrorCode (*monitordestroy)(void*),     &
      !        PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp     - iterative context obtained from KSPCreate()
      !   monitor - pointer to function (if this is PETSC_NULL, it
      !             turns off monitoring
      !   mctx    - [optional] context for private data for the monitor
      !             routine (use PETSC_NULL if no context is desired)
      !   monitordestroy - [optional] routine that frees monitor context
      !                    (may be PETSC_NULL)
      !
      ! Calling Sequence of monitor
      !
      ! call monitor (KSP ksp, int it, PetscReal rnorm, void *mctx)
      !
      !   ksp   - iterative context obtained from KSPCreate()
      !   it    - iteration number
      !   rnorm - (estimated) 2-norm of (preconditioned) residual
      !   mctx  - optional monitoring context, as set by KSPSetMonitor()
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetMonitor.html
      ! or PETSc users manual, pp.68

      ! i) default monitor (l_2 norm of the residual)
      ! call KSPSetMonitor(ksp,KSPDefaultMonitor, PETSC_NULL_OBJECT, &
      !                    PETSC_NULL_FUNCTION, PETScIerr)

      ! ii) user-defined monitor

      !call KSPSetMonitor(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
      !                   PETSC_NULL_FUNCTION, PETScIerr)


      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPSetMonitor")


!setflags!


!
!     ******************************************************************
!     *                                                                *
!     * Set the preconditioner. The PETSc libraries supports the       *
!     * following preconditioners:                                     *
!     *      PCJACOBI    - Jabobi                                      *
!     *      PCBJACOBI   - Block Jacobi                                *
!     *      PCSOR       - SOR                                         *
!     *      PCEISENSTAT - SOR with Eisenstat trick                    *
!     *      PCICC       - incomplete Cholesky                         *
!     *      PCILU       - incomplete LU                               *
!     *      PCASM       - additive Schwartz                           *
!     *      PCKSP       - linear solver                               *
!     *      PCCOMPOSITE - combination of preconditioners              *
!     *      PCLU        - LU                                          *
!     *      PCCHOLESKY  - Cholesky                                    *
!     *      PCNONE      - no preconditioning                          *
!     *      PCMG        - multigrid preconditioning                   *
!     *                                                                *
!     * For more information, see PETSc users manual, pp.70-77.        *
!     *                                                                *
!     ******************************************************************

      ! By extracting the KSP and PC contexts from the KSP context,
      ! we can then directly directly call any KSP and PC routines
      ! to set various options (PETSc users manual, pp.64).

      ! KSPGetPC - Returns a pointer to the preconditioner context set
      !            with KSPSetPC().
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPGetPC(KSP ksp,PC *pc,PetscErrorCode ierr)
      !
      ! Not Collective
      !
      ! Input Parameters
      !   ksp - iterative context obtained from KSPCreate()
      !
      ! Output Parameter
      !   pc - preconditioner context
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGetPC.html

      call KSPGetPC(ksp, pc, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScKsp", "Error in KSPGetPC")

      ! Set the preconditioning method.

      ! PCSetType - Builds PC for a particular preconditioner.
      !
      ! Synopsis
      !
      ! #include "petscpc.h" 
      ! call PCSetType(PC pc, PCType type, PetscErrorCode ierr)
      !
      ! Collective on PC
      !
      ! Input Parameter
      !   pc   - the preconditioner context.
      !   type - a known method 
      !
      ! see .../petsc/docs/manualpages/PC/PCSetType.html

!**********************************
!start of ASM code
!**********************************
      if (ASM)then    !pctype
     !try to use the Additive Schwartz Method ...!

      call PCSetType( pc, PCASM, PETScIerr)

!     
!     !Setup the local indexing schemes to for the ASM preconditioner
!
!      nPenaltyLocal = internalCommPenMatch(level,sps)%nPenalty
!
!      !allocate(localASMRange(0:nNodeslocal+nPenaltyLocal))
!      print *,'before allocate',PetscRank,((nNodeslocal+nPenaltyLocal)*nw)-1
!      allocate(localASMRange(0:(((nNodeslocal+nPenaltyLocal)*nw)-1)))
!      print *,'after allocate',PetscRank ,shape(localASMRange)     
!
!      !loop over the local domains and store the global nodes and the penalty nodes
!      counter = 0
!	
!      do nn=1,nDom
!	 call setPointers(nn,1,1)
!	 !globalNode => flowDoms(nn,1,1)%globalNode
!	 !print *,'globalshape',shape(globalNode),PetscRank
!
!	 do k = 1,kl
!   	    do j = 1,jl
!               do i = 1,il
!                  do n = 1,nw
!		     !localASMRange(counter)=globalNode(i,j,k)
!		     localASMRange(counter)=globalNode(i,j,k)*nw+n-1
!		     !print *,'counter',counter,i,j,k
!		     counter = counter +1
!		  enddo
!	       enddo
!            enddo
!         enddo
!      enddo

!!
!!     ****************************************************************
!!     *                                                              *
!!     * Copying the indexes of the 1 to 1 matching penalty data.     *
!!     *                                                              *
!!     ****************************************************************
!!
!      ! Set some variables to make the code more readable and test if
!      ! 1 to 1 matching penalty data must be copied.!
!
!      nNode = internalCommPenMatch(level,sps)%nPenalty
!
!      test11Copy: if(nNode > 0) then!
!
!        block1 => internalCommPenMatch(level,sps)%donorBlock
!        block2 => internalCommPenMatch(level,sps)%recvBlock
!        sub1   => internalCommPenMatch(level,sps)%donorSubface
!        sub2   => internalCommPenMatch(level,sps)%recvSubface
!        ind1   => internalCommPenMatch(level,sps)%donorInd
 !       ind2   => internalCommPenMatch(level,sps)%recvInd
!        indPen => internalCommPenMatch(level,sps)%donorIndPenalty

!        ! Loop over the nodes.!
!
!        do ii=1,nNode!
!
!          ! Store the indices a bit easier.
!
!          b1 = block1(ii); b2 = block2(ii)
!          s1 = sub1(ii);   s2 = sub2(ii)
!          i1 = ind1(ii,1); j1 = ind1(ii,2)
!          i2 = ind2(ii,1); j2 = ind2(ii,2)
!          ir = indPen(ii)
!
!
!          ! The global node numbering variable.
!	  do n=1,nw
!
!             !localASMRange(counter)=flowDoms(b1,level,sps)%penaltyData(s1)%globalNodeBuf(ir,i1,j1)
!	     localASMRange(counter)=flowDoms(b1,level,sps)%penaltyData(s1)%globalNodeBuf(ir,i1,j1)*nw+n-1
!	     !print *,'counter',counter,ir,i1,j1
!	     counter = counter +1
!           
!          enddo     
!
!        enddo 
!
!      endif test11Copy
!	
!      !length = (nNodeslocal+nPenaltyLocal)
!      length = (nNodeslocal+nPenaltyLocal)*nw
!      
!      !generate the IS index set for each subdomain
!      !call ISCreateGeneral(PETSC_COMM_WORLD,length,localASMRange,is, PETScIerr)
!      call ISCreateGeneral(PETSC_COMM_SELF,length,localASMRange,is, PETScIerr)
!	
!      deallocate(localASMRange)
!      
!      !Setup ASM Subdomains
!      Nsub = 1!_intType
!      call PCASMSetLocalSubdomains( pc, Nsub,is, PETScIerr)
!      
!      if( PETScIerr/=0 ) &
 !         call terminate("createPETScMat", &
 !                        "Error in PCASMSetLocalSubdomains")!
!!
!
 !     !Setup KSP context before calling local subdomain ksp contexts
 !     call KSPSetUp(ksp, PETScIerr);
 !    
 !     !get the subdomain contexts
 !     call PCASMGetSubKSP( pc, nlocal,  first, subksp, PETScIerr )
 !    
 !     !Setup the subdomain preconditioner
 !     call KSPGetPC( subksp, subpc, PETScIerr )!
!
!      !Set subdomain preconditioner type
 !     call PCSetType( subpc, PCILU, PETScIerr )!
!
!      !Set the matrix ordering type
!      !call PCFactorSetMatOrderingtype( subpc, MATORDERING_RCM, PETScIerr )
 !     call PCFactorSetMatOrderingtype( subpc,  MATORDERING_NATURAL, PETScIerr )
!
!      !Set ILU Parameters
!      call PCFactorSetLevels( subpc, 2 , PETScIerr)! // set 1 level of fill
!      call PCFactorSetFill( subpc, 2.5, PETScIerr )

 !     !Set sub KSP's as preconditioners only
 !     call KSPSetType(subksp, KSPPREONLY, PETScIerr)
               
!
  !    ! Convergence monitoring (default = silent mode).

      ! KSPSetMonitor - Sets an ADDITIONAL function to be called at
      !               every iteration to monitor the residual/error etc.
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetMonitor(KSP ksp, &
      !        PetscErrorCode (*monitor)(KSP,PetscInt,PetscReal,void*),&
      !        void *mctx,PetscErrorCode (*monitordestroy)(void*),     &
      !        PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameters
      !   ksp     - iterative context obtained from KSPCreate()
      !   monitor - pointer to function (if this is PETSC_NULL, it
      !             turns off monitoring
      !   mctx    - [optional] context for private data for the monitor
      !             routine (use PETSC_NULL if no context is desired)
      !   monitordestroy - [optional] routine that frees monitor context
      !                    (may be PETSC_NULL)
      !
      ! Calling Sequence of monitor
      !
      ! call monitor (KSP ksp, int it, PetscReal rnorm, void *mctx)
      !
      !   ksp   - iterative context obtained from KSPCreate()
      !   it    - iteration number
      !   rnorm - (estimated) 2-norm of (preconditioned) residual
      !   mctx  - optional monitoring context, as set by KSPSetMonitor()
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetMonitor.html
      ! or PETSc users manual, pp.68

      ! i) default monitor (l_2 norm of the residual)
      ! call KSPSetMonitor(ksp,KSPDefaultMonitor, PETSC_NULL_OBJECT, &
      !                    PETSC_NULL_FUNCTION, PETScIerr)

      ! i) default monitor (l_2 norm of the residual)
      ! call KSPSetMonitor(subksp,KSPDefaultMonitor, PETSC_NULL_OBJECT, &
      !                    PETSC_NULL_FUNCTION, PETScIerr)

    
     ! call KSPSetMonitor(subksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
      !                   PETSC_NULL_FUNCTION, PETScIerr)

      !call KSPSetMonitor(ksp,KSPDefaultMonitor, PETSC_NULL_OBJECT, &
      !                    PETSC_NULL_FUNCTION, PETScIerr)
 !     call KSPMonitorSet(subksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
 !                        PETSC_NULL_FUNCTION, PETScIerr)
 !     call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
 !                        PETSC_NULL_FUNCTION, PETScIerr)

      !call KSPMonitorSet(ksp,KSPMonitorDefault, PETSC_NULL_OBJECT, &
      !                    PETSC_NULL_FUNCTION, PETScIerr)




 !     if( PETScRank==0 ) then
 !       call KSPGetTolerances(subksp, rTol, aTol, dTol, mIts, PETScIerr)
 !       if( PETScIerr/=0 ) &
 !         call terminate("createPETScKsp", "Error in KSPGetTolerances")
 !       write(*,20) rTol, aTol, dTol
 !       write(*,21) mIts
  !    endif



      elseif(BJACOBI)then
!**************************
! block jacobi
!********************
      !print *,' entering block jacobi'
      !Set preconditioner type
      call PCSetType( pc, PCBJACOBI, PETScIerr)

      Nsub = 1!_intType
      length = nCellsLocal*nw
      call PCBJacobiSetLocalBlocks(pc,Nsub,length,PETScIerr)

      !Setup KSP context before calling local subdomain ksp contexts
      call KSPSetUp(ksp, PETScIerr)

      !Setup local subcontexts
      call PCBJacobiGetSubKSP(pc,nlocal,first,subksp,PETScIerr)

      !Setup Local ILU precondtioner
      call KSPGetPC( subksp, subpc, PETScIerr )
      call PCSetType( subpc, PCILU, PETScIerr )
 
      !Set the matrix ordering
      call PCFactorSetMatOrderingtype( subpc, MATORDERING_RCM, PETScIerr )
      !call PCFactorSetMatOrderingtype( subpc,  MATORDERING_NATURAL, PETScIerr )

      !Set ILU parameters
      call PCFactorSetLevels( subpc, 2 , PETScIerr)!  set 1 level of fill
      !call PCFactorSetFill( subpc, 3, PETScIerr )

      !Set local contexts to preconditioner's only
      call KSPSetType(subksp, KSPPREONLY, PETScIerr)

      ! set tolerances on subksp (is this really needed???)
      call KSPSetTolerances(subksp, 1.e-8, adjAbsTol, adjDivTol, &
                            adjMaxIter, PETScIerr)

      !Set the convergence monitors
      call KSPMonitorSet(subksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
                         PETSC_NULL_FUNCTION, PETScIerr)
     ! call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
     !                    PETSC_NULL_FUNCTION, PETScIerr)

      call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
                         PETSC_NULL_FUNCTION, PETScIerr)

     ! call KSPMonitorSet(ksp,KSPMonitorDefault, PETSC_NULL_OBJECT, &
     !                    PETSC_NULL_FUNCTION, PETScIerr)
!*********************
! Single Processor ILU
!*********************

      else
      print *,'inpcliu'
      !Set basic ILU (for single processor only!)

      !Set the convergence monitor
      call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
                         PETSC_NULL_FUNCTION, PETScIerr)
      
      ILULevels = 1_intType

      !Set preconditioner type
      call PCSetType(pc, PCILU, PETScIerr)
      !call PCSetType(pc, PCLU, PETScIerr)
   
      !Set ILU Levels
      call PCFactorSetLevels(pc, ILULevels, PETScIerr)
   
      !Set matrix ordering
      !!call PCFactorSetMatOrderingtype(pc,MATORDERING_RCM, PETScIerr)
      call PCFactorSetMatOrderingtype(pc,MATORDERING_ND, PETScIerr)

      endif !pctype

      !if( PETScIerr/=0 ) &
      !  call terminate("createPETScKsp", "Error in PCSetType")

      ! PCGetType - Gets the PC method type and name (as a string)
      !             from the PC context.
      ! Synopsis
      !
      ! #include "petscpc.h" 
      ! call PCGetType(PC pc,PCType meth, PetscErrorCode ierr)
      !
      ! Not Collective
      !
      ! Input Parameter
      !   pc   - the preconditioner context
      !
      ! Output Parameter
      !   meth - name of preconditioner 
      !
      ! see .../petsc/docs/manualpages/PC/PCGetType.html

      if( PETScRank==0 ) then
        call PCGetType(pc, pcType, PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("createPETScKsp", "Error in PCGetType")
        write(*,30) pcType
      endif

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
!
!     ******************************************************************
!
      subroutine MyKSPMonitor(myKsp, n, rnorm, dummy, ierr)
!
!     ******************************************************************
!     *                                                                *
!     * This is a user-defined routine for monitoring the KSP          *
!     * iterative solvers. Instead of outputing the L2-norm at every   *
!     * iteration (default PETSc monitor), it only does it every       *
!     * 'adjMonStep' iterations.                                       *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      implicit none
!
!     Subroutine arguments.
!
      ! myKsp - Iterative context
      ! n     - Iteration number
      ! rnorm - 2-norm (preconditioned) residual value
      ! dummy - Optional user-defined monitor context (unused here)
      ! ierr  - Return error code

      real(kind=realType), pointer, dimension(:,:) :: myKsp
      integer(kind=intType) :: n, dummy, ierr
      real(kind=realType)   :: rnorm
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! Write the residual norm to stdout every adjMonStep iterations.

      if( mod(n,adjMonStep)==0 ) then
        if( PETScRank==0 ) write(*,10) n, rnorm
      end if

      ierr = 0

      ! Output format.

   10 format(i4,1x,'KSP Residual norm',1x,e10.4)

#endif

      end subroutine MyKSPMonitor
