!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Allgatherv.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-02-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, &
                                 recvcounts, displs, recvtype, comm,    &
                                 error)
!
!      ******************************************************************
!      *                                                                *
!      * The subroutine MPI_Allgatherv calls MPI_ALLgather, possibly    *
!      * with an offset, which performs the same task as MPI_Allgatherv *
!      * for 1 processor.                                               *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments
!
       integer :: sendcount, sendtype, recvtype, comm, error

       integer, dimension(*)   :: recvcounts, displs
       character, dimension(*) :: sendbuf, recvbuf
!
!      Local variables.
!
       integer   :: start_index
!
!      Function definition
!
       integer :: determine_size_entity
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Call MPI_Allgather with a possible offset in the receive buffer.

       start_index = 1 + displs(1)*determine_size_entity(recvtype)

       call MPI_Allgather(sendbuf, sendcount, sendtype,        &
                          recvbuf(start_index), recvcounts(1), &
                          recvtype, comm, error)

      end subroutine MPI_Allgatherv
