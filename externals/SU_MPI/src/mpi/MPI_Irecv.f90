!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Irecv.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Irecv(buf, count, datatype, source, tag, comm, &
                            request, error)
!
!      ******************************************************************
!      *                                                                *
!      * The program should be such that no point to point              *
!      * communication takes place to the processor itself and thus     *
!      * MPI_Irecv is never called in sequential mode.                  *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments.
!
       integer :: count, datatype, source, tag, comm, request, error

       character, dimension(*) :: buf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call su_mpi_terminate("MPI_Irecv", &
                             "The routine MPI_Irecv should not be &
                             &called in sequential mode")

       end subroutine MPI_Irecv
