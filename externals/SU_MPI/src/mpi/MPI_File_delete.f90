!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_File_delete.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-28-2005                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_File_delete(filename, info, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_File_delete calls a C-routine, because the MPI-IO          *
!      * functionality cannot be simulated in Fortran.                  *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer :: info, error
       character(len=*) :: filename
!
!      Local variables.
!
       integer :: c_error, fileLen
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Call the C-function and set error accordingly.

       fileLen = len_trim(filename)
       call mpi_file_delete_c(filename, fileLen, c_error)
       if(c_error == 0) then
         error = mpi_success
       else
         error = mpi_error
       endif

       end subroutine MPI_File_delete
