!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Bcast.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-10-2002                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Bcast(buffer, count, datatype, root, comm, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_Bcast does not do anything in sequential mode.             *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       character(len=*) :: buffer
       integer          :: count, datatype, root, comm, error
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       error = mpi_success

       end subroutine MPI_Bcast
