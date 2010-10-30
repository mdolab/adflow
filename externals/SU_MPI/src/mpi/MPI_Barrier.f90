!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Barrier.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-11-2002                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Barrier(comm, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_Barrier does not do anything in sequential mode.           *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer :: comm, error
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       error = mpi_success

       end subroutine MPI_Barrier
