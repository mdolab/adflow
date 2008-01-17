/*
       ******************************************************************
       *                                                                *
       * File:          global_io_var.h                                 *
       * Author:        Edwin van der Weide                             *
       * Starting date: 01-13-2003                                      *
       * Last modified: 01-13-2003                                      *
       *                                                                *
       ******************************************************************
*/

#ifndef global_io_var_h
#define global_io_var_h

/*
       ******************************************************************
       *                                                                *
       * global_io_var.h contains some global variables needed to       *
       * simulate the MPI-IO functionalities in sequential mode.        *
       *                                                                *
       ******************************************************************
*/

#include <stdio.h>

/* The file pointers and their availabilities. */

#define NMAX_FILE_OPEN  50
extern FILE *fp_global[NMAX_FILE_OPEN];
extern int fp_global_used[NMAX_FILE_OPEN];

#endif /* global_io_var_h */
