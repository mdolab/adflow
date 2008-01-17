!
!      ******************************************************************
!      *                                                                *
!      * File:          MPI_File_open.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-13-2003                                      *
!      * Last modified: 02-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine MPI_File_open(comm, filename, amode, info, fh, error)
!
!      ******************************************************************
!      *                                                                *
!      * MPI_File_open calls a C-routine, because the MPI-IO            *
!      * functionality cannot be simulated in Fortran.                  *
!      *                                                                *
!      ******************************************************************
!
       use su_mpi
       implicit none
!
!      Subroutine arguments
!
       integer          :: comm, amode, info, fh, error
       character(len=*) :: filename
!
!      Local variables.
!
       integer :: c_error, c_amode, fileLen
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the open mode and set the correct value for c_amode.

       if(amode/mpi_mode_rdonly == 1) then
         c_amode = 1
       else if(amode/mpi_mode_rdwr == 1) then
         c_amode = 2
       else if(amode/mpi_mode_wronly == 1) then
         c_amode = 3
       else if(amode/mpi_mode_append == 1) then
         c_amode = 4
       else
         call su_mpi_terminate("MPI_File_open", &
                               "Unknown file access mode")
       endif

       ! Call the C-routine.

       fileLen = len_trim(filename)
       call mpi_file_open_c(filename, fileLen, c_amode, fh, c_error)

       if(c_error == 0) then
         error = mpi_success
       else if(c_error == 1) then
         call su_mpi_terminate("MPI_File_open", &
                                "Not enough file pointers. Change &
                                &nmax_file_open in global_io_var.h")
       else
         error = mpi_error
       endif

       end subroutine MPI_File_open
