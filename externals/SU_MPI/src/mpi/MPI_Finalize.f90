!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Finalize.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-10-2002                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Finalize(error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_Finalize does not do anything in sequential mode.          *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer :: error
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       error = mpi_success

       end subroutine MPI_Finalize
