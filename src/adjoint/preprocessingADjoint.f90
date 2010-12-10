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
subroutine preprocessingADjoint
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

  !     Local variables.
  !
  integer(kind=intType) :: ierr,level
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !

  ! Test the discrete adjoint solver assertions.
  level = 1_intType

  call assertionsADjoint(level)
  
  ! Allocate the memory for the global cell indexing of the
  ! computational mesh, later used to assemble the global
  ! adjoint system of equations.

  call allocMemADjoint(level)
     
  ! Determine the global cell and Node numbering.
  call setGlobalCellsAndNodes(level)

end subroutine preprocessingADjoint
