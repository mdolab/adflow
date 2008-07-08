subroutine SUmb_finalize()
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
!       *              Finalize Function for the mpi4py version          *
!       *                                                                *
!       * Because mpi4py is a module not an interpreter, the MPI needs   *
!       * to be initiallized in this library for it to run properly in   *
!       * parallel. This routine runs that initialization.               *
!       * Here parameter file is the name of the file with the input     *
!       * parameters.                                                    *
!       *                                                                *
!       ******************************************************************

use communication
!include "mpif.h"

integer rank, size, ierr

!   Clean up the MPI. 
call MPI_FINALIZE(ierr)

print *,'SUmb MPI FInalized'

end subroutine SUmb_finalize
