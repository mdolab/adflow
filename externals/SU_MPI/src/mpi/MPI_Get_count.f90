!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Get_count.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Get_count(status, datatype, count, error)
!
!      ******************************************************************
!      *                                                                *
!      * The program should be such that no point to point              *
!      * communication takes place to the processor itself and thus     *
!      * MPI_Get_count is never called in sequential mode.              *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments
!
       integer               :: datatype, count, error
       integer, dimension(*) :: status
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call su_mpi_terminate("MPI_Get_count", &
                             "The routine MPI_Get_count should not be &
                             &called in sequential mode")

       end subroutine MPI_Get_count
