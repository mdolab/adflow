!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Alltoall.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-20-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, &
                               recvcount, recvtype, comm, error)
!
!      ******************************************************************
!      *                                                                *
!      * The result of the subroutine MPI_Alltoall on 1 processor is a  *
!      * copy of the send buffer into the receive buffer.               *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer                 :: sendcount, sendtype
       integer                 :: recvcount, recvtype, comm, error
       character, dimension(*) :: sendbuf, recvbuf
!
!      Local variables
!
       integer :: i, nbytes_copy
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
       ! Determine the number of bytes to be copied and do so.

       nbytes_copy = recvcount*determine_size_entity(recvtype)
       do i=1,nbytes_copy
         recvbuf(i) = sendbuf(i)
       enddo

       ! Set error to mpi_success to indicate that everything
       ! went okay.

       error = mpi_success

       end subroutine MPI_Alltoall
