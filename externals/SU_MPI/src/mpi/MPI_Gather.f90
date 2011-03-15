!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Gather.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-14-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, &
                             recvcount, recvtype, root, comm, error)
!
!      ******************************************************************
!      *                                                                *
!      * In sequential mode the result of MPI_Gather is a copy of the   *
!      * send buffer into the receive buffer.                           *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer :: sendcount, sendtype, recvcount, recvtype 
       integer :: root, comm, error
       character, dimension(*) :: sendbuf, recvbuf
!
!      Local variables.
!
       integer :: i, nbytes_send, nbytes_recv, nbytes_copy
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
       ! Determine the sizes of the send and the receive buffer and
       ! take the minimum of the two.

       nbytes_send = sendcount*determine_size_entity(sendtype)
       nbytes_recv = recvcount*determine_size_entity(recvtype)
       nbytes_copy = min(nbytes_send, nbytes_recv)

       ! Copy the bytes from the send to the receive buffer.

       do i=1,nbytes_copy
         recvbuf(i) = sendbuf(i)
       enddo

       ! And set error to mpi_success to indicate that everything
       ! went okay.

       error = mpi_success

       end subroutine MPI_Gather
