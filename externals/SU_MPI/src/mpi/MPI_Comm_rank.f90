!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Comm_rank.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-10-2002                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Comm_rank(comm, rank, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_Comm_rank determines my processor ID in the group of the   *
!      * given communicator. In sequential mode this is simply 0.       *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer :: comm, rank, error
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       rank  = 0
       error = mpi_success

       end subroutine MPI_Comm_rank
