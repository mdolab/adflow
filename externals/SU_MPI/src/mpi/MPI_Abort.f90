!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Abort.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-10-2002                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Abort(comm, errorcode, error)
!
!      ******************************************************************
!      *                                                                *
!      * The subroutine MPI_Abort terminates the program.               *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer :: comm, errorcode, error
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       error = mpi_success
       stop

       end subroutine MPI_Abort
