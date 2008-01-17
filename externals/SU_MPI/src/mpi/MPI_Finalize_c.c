/*
       ******************************************************************
       *                                                                *
       * File:          MPI_Finalize_c.c                                *
       * Author:        Edwin van der Weide                             *
       * Starting date: 04-14-2004                                      *
       * Last modified: 02-21-2006                                      *
       *                                                                *
       ******************************************************************
*/

#include <stdio.h>

/*
       ******************************************************************
       *                                                                *
       * The function MPI_Finalize does not do anything in sequential   *
       * mode.                                                          *
       *                                                                *
       ******************************************************************
*/

int MPI_Finalize(void)
{
  return 0;
}
