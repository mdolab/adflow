!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Probe.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Probe(source, tag, comm, status, error)
!
!      ******************************************************************
!      *                                                                *
!      * The program should be such that no point to point              *
!      * communication takes place to the processor itself and thus     *
!      * MPI_Probe is never called in sequential mode.                  *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments.
!
       integer               :: source, tag, comm, error
       integer, dimension(*) :: status
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call su_mpi_terminate("MPI_Probe", &
                             "The routine MPIProbe should not be &
                             &called in sequential mode")

       end subroutine MPI_Probe
