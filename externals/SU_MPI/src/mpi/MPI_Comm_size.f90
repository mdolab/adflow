!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Comm_size.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-10-2002                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Comm_size(comm, size, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_Comm_size determines the number of processors in the group *
!      * of the given communicator. In sequential mode this is 1.       *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer :: comm, size, error
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       size  = 1
       error = mpi_success

       end subroutine MPI_Comm_size
