/*
       ******************************************************************
       *                                                                *
       * File:          MPI_File_delete_c.c                             *
       * Author:        Edwin van der Weide                             *
       * Starting date: 05-28-2005                                      *
       * Last modified: 02-21-2006                                      *
       *                                                                *
       ******************************************************************
*/

#include <stdio.h>
#include <stdlib.h>

/*
       ******************************************************************
       *                                                                *
       * Function which deletes the given file, if possible.            *
       *                                                                *
       ******************************************************************
*/

void MPI_File_delete_c(char *fileName, int *fileLen, int *error)
{
  int i;
  char command[1024];

  /* Add the terminating character to fileName at the correct place. */

  i = *fileLen;
  fileName[i] = '\0';

  /* Create the command for deleting the file and do the system call. */
  /* Set error to 0 to indicate that everything went okay.            */

  sprintf(command, "rm -f %s", fileName);
  system(command);
  error = 0;
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void MPI_FILE_DELETE_C(char *fileName, int *fileLen, int *error)
{
  MPI_File_delete_c(fileName, fileLen, error);
}

void mpi_file_delete_c(char *fileName, int *fileLen, int *error)
{
  MPI_File_delete_c(fileName, fileLen, error);
}

void mpi_file_delete_c_(char *fileName, int *fileLen, int *error)
{
  MPI_File_delete_c(fileName, fileLen, error);
}

void mpi_file_delete_c__(char *fileName, int *fileLen, int *error)
{
  MPI_File_delete_c(fileName, fileLen, error);
}
