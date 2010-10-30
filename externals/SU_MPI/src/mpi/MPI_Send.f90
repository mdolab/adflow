!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Send.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Send(buf, count, datatype, dest, tag, comm, error)
!
!      ******************************************************************
!      *                                                                *
!      * The program should be such that no point to point              *
!      * communication takes place to the processor itself and thus     *
!      * MPI_Send is never called in sequential mode.                   *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments.
!
       integer :: count, datatype, dest, tag, comm, error

       character, dimension(*) :: buf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call su_mpi_terminate("MPI_Send", &
                             "The routine MPI_Send should not be &
                             &called in sequential mode")

       end subroutine MPI_Send
