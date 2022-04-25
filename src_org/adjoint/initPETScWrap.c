/*
       ******************************************************************
       *                                                                *
       * File:          initPETScWrap.c                                 *
       * Author:        Andre C. Marta                                  *
       * Starting date: 02-27-2007                                      *
       * Last modified: 02-27-2007                                      *
       *                                                                *
       ******************************************************************
*/

#ifndef USE_NO_PETSC

#include "petsc.h"

/*
       ******************************************************************
       *                                                                *
       * The function initPETSc call the C interface that initializes   *
       * the PETSc library. This was necessary because a call to the    *
       * fortran PETSc interface causes a segmentation fault for some   *
       * unknown reason.                                                *
       *                                                                *
       * The original initialization PETSc C function is coded in       *
       * <petsc dir>/src/sys/objects/pinit.c line 396--                 *
       *                                                                *
       * Fortran usage of this wrapped routine:                         *
       * >> call initPETScWrap()                                        *
       *                                                                *
       ******************************************************************
*/

void initPETScWrap(void)
{
  int ierr;
  int argc;

  argc=0;
    ierr = PetscInitialize(&argc, (char ***)0, (char *)0, (char *)0);

}

/*
       ******************************************************************
       *                                                                *
       * Dummy functions to allow the calling from FORTRAN.             *
       *                                                                *
       ******************************************************************
*/

void INITPETSCWRAP(void)
{
  initPETScWrap();
}

void initpetscwrap_(void)
{
  initPETScWrap();
}

void initpetscwrap__(void)
{
  initPETScWrap();
}

#endif
