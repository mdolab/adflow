/*
       ******************************************************************
       *                                                                *
       * File:          my_time.c                                       *
       * Author:        Edwin van der Weide                             *
       * Starting date: 11-26-2003                                      *
       * Last modified: 07-10-2005                                      *
       *                                                                *
       ******************************************************************
*/

#include <stdio.h>
#include <time.h>

/*
       ******************************************************************
       *                                                                *
       * Function to compute the time relative to a time in the past.   *
       *                                                                *
       ******************************************************************
*/

double my_time(void)
{
  clock_t elapsedTime = clock();
  return ((double) elapsedTime) / ((double) CLOCKS_PER_SEC);
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

double MY_TIME(){return my_time();}
double my_time_(){return my_time();}
double my_time__(){return my_time();}
