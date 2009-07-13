!
!     ******************************************************************
!     *                                                                *
!     * File:          preprocessingADjoint.f90                        *
!     * Author:        Andre C. Marta , C.A.(Sandy) Mader              *
!     *                Seongim Choi
!     * Starting date: 07-20-2006                                      *
!     * Last modified: 01-18-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine preprocessingADjoint(level)
!
!     ******************************************************************
!     *                                                                *
!     * Perform the preprocessing tasks for the adjoint solver, such   *
!     * as assertions testing, memory allocation and global indexing.  *
!     *                                                                *
!     ******************************************************************
!
      use communication
      use precision
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level
!
!     Local variables.
!
      integer :: ierr
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      !print *,' in preprocessing adjoint'
      ! Test the discrete adjoint solver assertions.

      call assertionsADjoint(level)
      !print *,'assertions called'
      ! Allocate the memory for the global cell indexing of the
      ! computational mesh, later used to assemble the global
      ! adjoint system of equations.

      call allocMemADjoint(level)
      !print *,'memory allocated'
      ! Determine the global cell and Node numbering.

      call setGlobalCellsAndNodes(level)
      print *,'global node indices set'
      ! Synchronize the processors.

      call mpi_barrier(SUmb_comm_world, ierr)

      end subroutine preprocessingADjoint
