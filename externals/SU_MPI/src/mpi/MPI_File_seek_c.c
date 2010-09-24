/*
       ******************************************************************
       *                                                                *
       * File:          MPI_File_seek_c.c                               *
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
       * Function that sets the file pointer to the new position.       *
       *                                                                *
       ******************************************************************
*/

void MPI_File_seek_c(int *fh, long long *offset, int *whence, int *error)
{
  /* Determine the case and sets the file pointer accordingly. */

  switch( *whence )
  {
    case 1:

      /* Position the file pointer relative from the beginning. */

      *error = fseek(fp_global[*fh], *offset, SEEK_SET);
      break;

    case 2:

      /* Position the file pointer relative from the current position. */

      *error = fseek(fp_global[*fh], *offset, SEEK_CUR);
      break;

    case 3:

      /* Position the file pointer relative from the end. */

      *error = fseek(fp_global[*fh], *offset, SEEK_END);
      break;
  }
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void MPI_FILE_SEEK_C(int *fh, long long *offset, int *whence, int *error)
{
  MPI_File_seek_c(fh, offset, whence, error);
}

void mpi_file_seek_c(int *fh, long long *offset, int *whence, int *error)
{
  MPI_File_seek_c(fh, offset, whence, error);
}

void mpi_file_seek_c_(int *fh, long long *offset, int *whence, int *error)
{
  MPI_File_seek_c(fh, offset, whence, error);
}

void mpi_file_seek_c__(int *fh, long long *offset, int *whence, int *error)
{
  MPI_File_seek_c(fh, offset, whence, error);
}
