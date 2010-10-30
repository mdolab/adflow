/*
       ******************************************************************
       *                                                                *
       * File:          f77flush.c                                      *
       * Author:        Andre C. Marta                                  *
       * Starting date: 10-06-2006                                      *
       * Last modified: 10-06-2006                                      *
       *                                                                *
       ******************************************************************
*/

#include <stdio.h>

/*
       ******************************************************************
       *                                                                *
       * The function f77flush flushes the Fortran output buffer. This  *
       * function only exists in C, therefore a call to C from Fortran  *
       * is required:                                                   *
       *                                                                *
       * >> call f77flush()                                             *
       *                                                                *
       ******************************************************************
*/

void f77flush(void)
{
fflush(stdout);
fflush(stderr);
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void F77FLUSH(void)
{
  f77flush();
}

void f77flush_(void)
{
  f77flush();
}

void f77flush__(void)
{
  f77flush();
}
