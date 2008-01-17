!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Allreduce.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-13-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, &
                                comm, error)
!
!      ******************************************************************
!      *                                                                *
!      * The subroutine MPI_Allreduce calls MPI_Reduce, which performs  *
!      * the same task as MPI_Allreduce for 1 processor.                *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments
!
       integer                 :: count, datatype, op, comm, error
       character, dimension(*) :: sendbuf, recvbuf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call MPI_Reduce(sendbuf, recvbuf, count, datatype, op, &
                       0, comm, error)

       end subroutine MPI_Allreduce
