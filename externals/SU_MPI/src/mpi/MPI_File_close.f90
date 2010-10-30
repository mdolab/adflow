!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_File_close.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-13-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_File_close(fh, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_File_close calls a C-routine, because the MPI-IO           *
!      * functionality cannot be simulated in Fortran.                  *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer :: fh, error
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Call the C-function and set error to mpi_success

       call mpi_file_close_c(fh)
       error = mpi_success

       end subroutine MPI_File_close
