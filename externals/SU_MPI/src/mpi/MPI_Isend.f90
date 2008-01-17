!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Isend.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Isend(buf, count, datatype, dest, tag, comm, &
                            request, error)
!
!      ******************************************************************
!      *                                                                *
!      * The program should be such that no point to point              *
!      * communication takes place to the processor itself and thus     *
!      * MPI_Isend is never called in sequential mode.                  *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments.
!
       integer :: count, datatype, dest, tag, comm, request, error

       character, dimension(*) :: buf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call su_mpi_terminate("MPI_Isend", &
                             "The routine MPI_Isend should not be &
                             &called in sequential mode")

       end subroutine MPI_Isend
