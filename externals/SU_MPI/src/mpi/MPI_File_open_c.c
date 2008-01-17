/*
       ******************************************************************
       *                                                                *
       * File:          MPI_File_open_c.c                               *
       * Author:        Edwin van der Weide                             *
       * Starting date: 01-13-2003                                      *
       * Last modified: 04-21-2006                                      *
       *                                                                *
       ******************************************************************
*/

#include <string.h>
#include "global_io_var.h"

/*
       ******************************************************************
       *                                                                *
       * Definition of the global variables needed for the file IO.     *
       *                                                                *
       ******************************************************************
*/

FILE *fp_global[NMAX_FILE_OPEN];
int fp_global_used[NMAX_FILE_OPEN];

/*
       ******************************************************************
       *                                                                *
       * Function which opens the file for the given mode.              *
       *                                                                *
       ******************************************************************
*/

void MPI_File_open_c(char *fileName, int *fileLen, int *aMode, int *fh,
                     int *error)
{
  int i;
  static int firstCall = 1;

  /* If this is the first call to MPI_File_open_c initialize         */
  /* fp_global_used to 0. It might be that globals are automatically */
  /* initialized to 0, so this loop might not be necessary. However  */
  /* I am not entirely sure and therefore fp_global_used is          */
  /* explicitly initialized here.                                    */

  if( firstCall )
  {
    for(i=0; i<NMAX_FILE_OPEN; i++) fp_global_used[i] = 0;
    firstCall = 0;
  }

  /* Add the terminating character to fileName at the correct place. */

  i = *fileLen;
  fileName[i] = '\0';

  /* Find the first available file pointer. */

  for(i=0; i<NMAX_FILE_OPEN; i++)
    if( !fp_global_used[i] ) break;

  /* Check if a file pointer can be used. If not set error to 1 and */
  /* return. An error value of 1 will cause the Fortran program to  */
  /* take care of the error; this is a bit more convenient.         */

  if(i == NMAX_FILE_OPEN)
  {
    *error = 1;
    return;
  }

  /* Try to open the file. */

  switch( *aMode )
  {
    case 1:

      /* Open for read only. */

      fp_global[i] = fopen(fileName, "r");
      break;

    case 2:

      /* Open for reading and writing. */

      fp_global[i] = fopen(fileName, "r+");
      break;

    case 3:

      /* Open for write only. */

      fp_global[i] = fopen(fileName, "w");
      break;

    case 4:

      /* Open for appending. */

      fp_global[i] = fopen(fileName, "a");
      break;
  }

  /* Check if it went okay. If not set error to 2 and return. */

  if( !fp_global[i] )
  {
    *error = 2;
    return;
  }

  /* Set the file handler fh to i, fp_global_used[i] to true and */
  /* error to 0 to indicate that everything went okay.           */

  *fh               = i;
  fp_global_used[i] = 1;
  *error            = 0;
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void MPI_FILE_OPEN_C(char *fileName, int *fileLen, int *aMode, int *fh,
                     int *error)
{
  MPI_File_open_c(fileName, fileLen, aMode, fh, error);
}

void mpi_file_open_c(char *fileName, int *fileLen, int *aMode, int *fh,
                     int *error)
{
  MPI_File_open_c(fileName, fileLen, aMode, fh, error);
}

void mpi_file_open_c_(char *fileName, int *fileLen, int *aMode, int *fh,
                      int *error)
{
  MPI_File_open_c(fileName, fileLen, aMode, fh, error);
}

void mpi_file_open_c__(char *fileName, int *fileLen, int *aMode, int *fh,
                       int *error)
{
  MPI_File_open_c(fileName, fileLen, aMode, fh, error);
}
