!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_Waitall.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_Waitall(count, array_of_requests, &
                              array_of_statuses, error)
!
!      ******************************************************************
!      *                                                                *
!      * The program should be such that no point to point              *
!      * communication takes place to the processor itself and thus     *
!      * MPI_Waitall should only be called in sequential mode with      *
!      * count = 0.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer               :: count, error
       integer, dimension(*) :: array_of_requests, array_of_statuses
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       if(count == 0) then

         ! Simply set error to mpi_success

         error = mpi_success

       else

         call su_mpi_terminate("MPI_Waitall", &
                               "The routine MPI_Waitall should not be &
                               &called in sequential mode with &
                               &non-zero requests to complete.")
       endif

       end subroutine MPI_Waitall
