/*
       ******************************************************************
       *                                                                *
       * File:          myIsNaNC.c                                      *
       * Author:        Edwin van der Weide                             *
       * Starting date: 10-21-2007                                      *
       * Last modified: 10-21-2007                                      *
       *                                                                *
       ******************************************************************
*/

#include "SUmb_c_types.h"
#include <math.h>

/*
       ******************************************************************
       *                                                                *
       * The function myIsNanC determines whether or not the the given  *
       * floating point value is a NaN or infinity. This function is    *
       * typically used when Fortran compilers are used, which do not   *
       * have a isnan function.                                         *
       *                                                                *
       ******************************************************************
*/

void myIsNanC(SUmb_floatT *val, SUmb_intT *res)
{
  if(isnan(*val) || isinf(*val)) *res = 1;
  else                           *res = 0;
}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void MYISNANC(SUmb_floatT *val, SUmb_intT *res) {myIsNanC(val, res);}
void myisnanc_(SUmb_floatT *val, SUmb_intT *res) {myIsNanC(val, res);}
void myisnanc__(SUmb_floatT *val, SUmb_intT *res) {myIsNanC(val, res);}
