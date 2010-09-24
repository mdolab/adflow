/*
       ******************************************************************
       *                                                                *
       * File:          SUmb_c_types.h                                  *
       * Author:        Edwin van der Weide                             *
       * Starting date: 02-11-2003                                      *
       * Last modified: 06-12-2005                                      *
       *                                                                *
       ******************************************************************
*/

#ifndef SUmb_c_types_h
#define SUmb_c_types_h

/*
       ******************************************************************
       *                                                                *
       * SUmb_c_types.h defines the integer and floating point types to *
       * be used for the C source code, depending on the compiler       *
       * options. In this way it is compatible with the kind type used  *
       * in the Fortran sources.                                        *
       *                                                                *
       ******************************************************************
*/

#include <stdio.h>

/****************************/
/* Integer type definition. */
/****************************/

#ifdef USE_LONG_INT

  /* 8 byte integers must be used for the default integer type. */

  typedef long long SUmb_intT;

#else

  /* Default, 4 byte, integers must be used for the default integer type. */

  typedef int SUmb_intT;

#endif

/***********************************/
/* Floating point type definition. */
/***********************************/

#ifdef USE_SINGLE_PRECISION

  /* Single precision floating point type must be used. */

  typedef float  SUmb_floatT;

#elif USE_QUADRUPLE_PRECISION

  /* Quadruple precision floating point type must be used. */
  /* This may not be fully compatible for 32 bit machines. */

  typedef long double SUmb_floatT;

#else

  /* Double precision floating point type must be used. */

  typedef double SUmb_floatT;

#endif

#endif /* SUmb_c_types_h */
