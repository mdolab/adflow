!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Recv.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Recv(buf, count, datatype, source, tag, comm, &
                           status, error)
!
!      ******************************************************************
!      *                                                                *
!      * The program should be such that no point to point              *
!      * communication takes place to the processor itself and thus     *
!      * MPI_Recv is never called in sequential mode.                   *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments.
!
       integer               :: count, datatype, source, tag, comm, error
       integer, dimension(*) :: status

       character, dimension(*) :: buf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call su_mpi_terminate("MPI_Recv", &
                             "The routine MPI_Recv should not be &
                             &called in sequential mode")

       end subroutine MPI_Recv
