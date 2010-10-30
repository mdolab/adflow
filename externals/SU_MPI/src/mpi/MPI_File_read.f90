!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_File_read.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-13-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_File_read(fh, buf, count, datatype, status, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_File_read calls a C-routine, because the MPI-IO            *
!      * functionality cannot be simulated in Fortran.                  *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer                 :: fh, count, datatype, error
       integer, dimension(*)   :: status
       character, dimension(*) :: buf
!
!      Local variables.
!
       integer :: size_entity, c_error
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
       ! Determine the size of each individual entity.

       size_entity = determine_size_entity(datatype)

       ! Call the C-routine.

       call mpi_file_read_c(fh, buf, count, size_entity, c_error)

       if(c_error == 0) then
         error = mpi_success
       else
         error = mpi_error
       endif

       end subroutine MPI_File_read
