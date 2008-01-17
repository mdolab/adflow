!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Waitany.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Waitany(count, array_of_requests, index, status, &
                              error)
!
!      ******************************************************************
!      *                                                                *
!      * The program should be such that no point to point              *
!      * communication takes place to the processor itself and thus     *
!      * MPI_Waitany is never called in sequential mode.                *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments
!
       integer               :: count, index, error
       integer, dimension(*) :: array_of_requests, status
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call su_mpi_terminate("MPI_Waitany", &
                             "The routine MPI_Waitany should not be &
                             &called in sequential mode")

       end subroutine MPI_Waitany
