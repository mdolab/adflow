/*
       ******************************************************************
       *                                                                *
       * File:          MPI_File_close_c.c                              *
       * Author:        Edwin van der Weide                             *
       * Starting date: 01-13-2003                                      *
       * Last modified: 03-06-2006                                      *
       *                                                                *
       ******************************************************************
*/

#include "global_io_var.h"

/*
       ******************************************************************
       *                                                                *
       * Function which closes the given file.                          *
       *                                                                *
       ******************************************************************
*/

void MPI_File_close_c(int *fh)
{
  /* Close the file and resets the corresponding value of fp_global_used. */

  fclose(fp_global[*fh]);
  fp_global_used[*fh] = 0;
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void MPI_FILE_CLOSE_C(int *fh)
{
  MPI_File_close_c(fh);
}

void mpi_file_close_c(int *fh)
{
  MPI_File_close_c(fh);
}

void mpi_file_close_c_(int *fh)
{
  MPI_File_close_c(fh);
}

void mpi_file_close_c__(int *fh)
{
  MPI_File_close_c(fh);
}
