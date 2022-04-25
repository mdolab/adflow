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

#include "ADflow_c_types.h"
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

void myIsNanC(ADflow_floatT *val, ADflow_intT *res)
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

void MYISNANC(ADflow_floatT *val, ADflow_intT *res) {myIsNanC(val, res);}
void myisnanc_(ADflow_floatT *val, ADflow_intT *res) {myIsNanC(val, res);}
void myisnanc__(ADflow_floatT *val, ADflow_intT *res) {myIsNanC(val, res);}
