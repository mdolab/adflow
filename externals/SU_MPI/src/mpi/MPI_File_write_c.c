/*
       ******************************************************************
       *                                                                *
       * File:          MPI_File_write_c.c                              *
       * Author:        Edwin van der Weide                             *
       * Starting date: 03-20-2005                                      *
       * Last modified: 02-21-2006                                      *
       *                                                                *
       ******************************************************************
*/

#include "global_io_var.h"

/*
       ******************************************************************
       *                                                                *
       * Function which writes a number of bytes to the given file.     *
       *                                                                *
       ******************************************************************
*/

void MPI_File_write_c(int *fh, void *buf, int *count, int *size,
                      int *error)
{
  size_t nItemsWritten;

  /* Write the items to the stream. */

  nItemsWritten = fwrite(buf, *size, *count, fp_global[*fh]);

  /* Set error, depending on the result of the reading. */

  if(nItemsWritten == *count) *error = 0;
  else                        *error = 2;
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void MPI_FILE_WRITE_C(int *fh, void *buf, int *count, int *size,
                      int *error)
{
  MPI_File_write_c(fh, buf, count, size, error);
}

void mpi_file_write_c_(int *fh, void *buf, int *count, int *size,
                       int *error)
{
  MPI_File_write_c(fh, buf, count, size, error);
}

void mpi_file_write_c__(int *fh, void *buf, int *count, int *size,
                        int *error)
{
  MPI_File_write_c(fh, buf, count, size, error);
}
