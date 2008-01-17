!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Comm_split.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-20-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Comm_split(comm, color, key, newcomm, error)
!
!      ******************************************************************
!      *                                                                *
!      * In sequential mode MPI_Comm_split only needs to copy the       *
!      * communicator.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer :: comm, color, key, newcomm, error
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       newcomm = comm
       error = mpi_success

       end subroutine MPI_Comm_split
