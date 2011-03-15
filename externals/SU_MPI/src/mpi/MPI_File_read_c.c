/*
       ******************************************************************
       *                                                                *
       * File:          MPI_File_read_c.c                               *
       * Author:        Edwin van der Weide                             *
       * Starting date: 01-13-2003                                      *
       * Last modified: 02-21-2006                                      *
       *                                                                *
       ******************************************************************
*/

#include "global_io_var.h"

/*
       ******************************************************************
       *                                                                *
       * Function which reads a number of bytes from the given file.    *
       *                                                                *
       ******************************************************************
*/

void MPI_File_read_c(int *fh, void *buf, int *count, int *size,
                     int *error)
{
  size_t nItemsRead;

  /* Read the items from the stream. */

  nItemsRead = fread(buf, *size, *count, fp_global[*fh]);

  /* Set error, depending on the result of the reading. */

  if(nItemsRead == *count) *error = 0;
  else                     *error = 2;
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void MPI_FILE_READ_C(int *fh, void *buf, int *count, int *size, int *error)
{
  MPI_File_read_c(fh, buf, count, size, error);
}

void mpi_file_read_c(int *fh, void *buf, int *count, int *size, int *error)
{
  MPI_File_read_c(fh, buf, count, size, error);
}

void mpi_file_read_c_(int *fh, void *buf, int *count, int *size, int *error)
{
  MPI_File_read_c(fh, buf, count, size, error);
}

void mpi_file_read_c__(int *fh, void *buf, int *count, int *size, int *error)
{
  MPI_File_read_c(fh, buf, count, size, error);
}
