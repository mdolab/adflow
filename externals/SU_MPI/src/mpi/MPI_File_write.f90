!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_File_write.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-20-2005                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_File_write(fh, buf, count, datatype, status, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_File_write calls a C-routine, because the MPI-IO           *
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

       ! Call the C-routine and set the error.

       call mpi_file_write_c(fh, buf, count, size_entity, c_error)

       if(c_error == 0) then
         error = mpi_success
       else
         error = mpi_error
       endif

       end subroutine MPI_File_write
