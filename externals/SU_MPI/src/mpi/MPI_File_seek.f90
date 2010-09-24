!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_File_seek.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-13-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_File_seek(fh, offset, whence, error)
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
       integer                       :: fh, whence, error
       integer(kind=mpi_offset_kind) :: offset
!
!      Local variables.
!
       integer         :: c_whence, c_error
       integer(kind=8) :: c_offset
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the value of c_whence depending on whence.

       select case (whence)
         case (mpi_seek_set)
           c_whence = 1
         case (mpi_seek_cur)
           c_whence = 2
         case (mpi_seek_end)
           c_whence = 3
       end select

       ! Copy offset in c_offset to be sure that an integer8
       ! (long long) is sent to the C function.

       c_offset = offset

       ! Call the C-function and set error.

       call mpi_file_seek_c(fh, c_offset, c_whence, c_error)
       if(c_error == 0) then
         error = mpi_success
       else 
         error = mpi_error
       endif

       end subroutine MPI_File_seek
