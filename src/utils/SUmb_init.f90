subroutine SUmb_init()
!
!       ******************************************************************
!       *                                                                *
!       * File:          SUmb_init.f90                                   *
!       * Author:        C.A.(Sandy) Mader                               *
!       * Starting date: 10-15-2007                                      *
!       * Last modified: 07-02-2008                                      *
!       *                                                                *
!       ******************************************************************

!
!       ******************************************************************
!       *                                                                *
!       *              Init Function for the mpi4py version              *
!       *                                                                *
!       * Because mpi4py is a module not an interpreter, the MPI needs   *
!       * to be initiallized in this library for it to run properly in   *
!       * parallel. This routine runs that initialization.               *
!       * Here parameter file is the name of the file with the input     *
!       * parameters.                                                    *
!       *                                                                *
!       ******************************************************************

use communication
implicit none
!include "mpif.h"

integer rank, size, ierr

!   Enroll the program in MPI; for a sequential code this is a dummy. 
!call MPI_INIT(ierr)

call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

if( rank ==0)then
   print *,'SUmb MPI Initialization Parameters'
   print *,'MPI_COMM_WORLD = ',MPI_COMM_WORLD
   print *,'Size = ',size
endif
print *,'Rank = ',rank
SUmb_comm_world = mpi_comm_world
SUmb_comm_self = mpi_comm_self
end subroutine SUmb_init
