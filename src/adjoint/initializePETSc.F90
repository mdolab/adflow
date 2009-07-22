!
!     ******************************************************************
!     *                                                                *
!     * File:          initializePETSc.F90                             *
!     * Author:        C.A.(Sandy) Mader, Andre C. Marta               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine initializePETSc
!
!     ******************************************************************
!     *                                                                *
!     * Initialize PETSc by calling the appropriate routine provided   *
!     * provided in the PETSc library - PetscInitialize. This          *
!     * automatically calls MPI_Init() if MPI has not been previously  *
!     * initialized. It also initializes some auxiliar variables       *
!     * declared in module ADjointPETSc.                               *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use communication
      implicit none
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! The PETSc initialization, which might call mpi_init, has to be
      ! placed right after where the call to mpi_init in SUmb is, that
      ! is to say, it has moved to SUmb.c

      ! PetscInitialize - Initializes the PETSc database and MPI.
      !   PetscInitialize() calls MPI_Init() if that has yet to be
      !   called, so this routine should always be called near the
      !   beginning of your program -- usually the very first line!
      !
      ! Synopsis
      !
      ! #include "petsc.h"   
      ! call PetscInitialize(int *argc,char ***args,const char file[], &
      !                      const char help[],PetscErrorCode ierr)
      !
      ! Collective on MPI_COMM_WORLD or PETSC_COMM_WORLD
      !   if it has been set
      !
      ! Input Parameters
      !   argc - count of number of command line arguments
      !   args - the command line arguments
      !   file - [optional] PETSc database file, defaults to
      !           ~username/.petscrc (use PETSC_NULL for default)
      !   help - [optional] Help message to print, use PETSC_NULL
      !          for no message
      !
      ! If you wish PETSc to run on a subcommunicator of MPI_COMM_WORLD,
      !   create that communicator first and assign it to
      !   PETSC_COMM_WORLD BEFORE calling PetscInitialize() 
      !
      ! Notes
      ! If for some reason you must call MPI_Init() separately, call it
      !   before PetscInitialize().
      !
      ! Fortran Version
      ! In Fortran this routine has the format
      !
      ! call PetscInitialize(file,ierr)
      !
      !   ierr - error return code
      !   file - [optional] PETSc database file name, defaults to
      !        ~username/.petscrc (use PETSC_NULL_CHARACTER for default)
      ! Important Fortran Note
      ! In Fortran, you MUST use PETSC_NULL_CHARACTER to indicate a null
      !   character string; you CANNOT just use PETSC_NULL as in the C
      !   version. See the users manual for details.
      !
      ! see .../petsc/docs/manualpages/Sys/PetscInitialize.html
      ! or PETSc users manual, pp.17

      ! call PetscInitialize(PETSC_NULL_CHARACTER, PETScIerr)

      ! Note: the above routine causes a segmentation fault
      !       for unknown reasons. A work around was to call
      !       the C-version, instead of the Fortran, which
      !       as been wrapped in initPETScWrap.c

      call initPETScWrap()

      !setup debugger
      !call PetscAttachDebugger(PETScIerr)
!      call PetscMallocDebug(PETSC_TRUE, PETScIerr)
!      call PetscOptionsSetValue('-malloc_debug','' ,PETScIerr)
!      print *, "PetscOptionsSetValue PETScIerr =", PETScIerr

      ! Set the PETSc communicator to the SUmb communicator.

      PETSC_COMM_WORLD = SUmb_comm_world
      PETSC_COMM_SELF = SUmb_comm_self

      ! Determine the communicator size and the processor rank.

      call MPI_Comm_size(PETSC_COMM_WORLD, PETScSize, PETScIerr)
      call MPI_Comm_rank(PETSC_COMM_WORLD, PETScRank, PETScIerr)

      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Initializing PETSc..."

      ! Initialize some auxiliar variables in module ADjointPETSc.

      PETScNegOne = -1.0_realType
      PETScZero   =  0.0_realType
      PETScOne    =  1.0_realType

      ! Allocate memory for the convergence residual history.

      if (.not. allocated(adjResHist))then
         allocate(adjResHist(adjMaxIter))
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output format.

   10 format(a)

#endif

      end subroutine initializePETSc
