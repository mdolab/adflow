!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Allgather.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, &
                                recvcount, recvtype, comm, error)
!
!      ******************************************************************
!      *                                                                *
!      * The subroutine MPI_Allgather calls MPI_Gather, which performs  *
!      * the same task as MPI_Allgather for 1 processor.                *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments
!
       integer :: sendcount, sendtype, recvcount, recvtype, comm, error
       character, dimension(*) :: sendbuf, recvbuf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, &
                       recvcount, recvtype, 0, comm, error)

       end subroutine MPI_Allgather
