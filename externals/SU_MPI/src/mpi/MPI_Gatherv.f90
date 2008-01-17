!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Gatherv.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2004                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, &
                              recvcounts, displs, recvtype, root,    &
                              comm, error)
!
!      ******************************************************************
!      *                                                                *
!      * The result of MPI_Gatherv in sequential mode is identical to   *
!      * MPI_Gather with a possible offset.                             *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments
!
       integer :: sendcount, sendtype, recvtype, root, comm, error

       integer, dimension(*)   :: recvcounts, displs
       character, dimension(*) :: sendbuf, recvbuf
!
!      Local variables.
!
       integer :: start_index
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
       ! Call MPI_Gather with a possible offset in the receive buffer.

       start_index = 1 + displs(1)*determine_size_entity(recvtype)

       call MPI_Gather(sendbuf, sendcount, sendtype,        &
                       recvbuf(start_index), recvcounts(1), &
                       recvtype, root, comm, error)

       end subroutine MPI_Gatherv
