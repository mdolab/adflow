/*
       ******************************************************************
       *                                                                *
       * File:          SUmb.c                                          *
       * Author:        Edwin van der Weide                             *
       * Starting date: 04-14-2004                                      *
       * Last modified: 07-07-2005                                      *
       *                                                                *
       ******************************************************************
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "SUmb_c_types.h"

/*
       ******************************************************************
       *                                                                *
       *                  Function prototypes.                          *
       *                                                                *
       ******************************************************************
*/

int MPI_Init(int *argc, char ***argv);
int MPI_Finalize(void);

#ifdef FORTRAN_CAPITALS

  void SUMB_FORTRAN(const char *exec_name,  const int *len_exec,
		    const char *param_name, const int *len_param,
		    const int *nargs,       const int  *size_int,
                    int len1,               int len2);

#elif FORTRAN_DOUBLE_UNDERSCORE

  void sumb_fortran__(const char *exec_name,  const int *len_exec,
		      const char *param_name, const int *len_param,
		      const int *nargs,       const int  *size_int,
                      int len1,               int len2);

#elif FORTRAN_NO_UNDERSCORE

  void sumb_fortran(const char *exec_name,  const int *len_exec,
                    const char *param_name, const int *len_param,
                    const int *nargs,       const int  *size_int,
                    int len1,               int len2);

#else

  void sumb_fortran_(const char *exec_name,  const int *len_exec,
		     const char *param_name, const int *len_param,
		     const int *nargs,       const int  *size_int,
                     int len1,               int len2);

#endif

/*
       ******************************************************************
       *                                                                *
       *              Main program of the SUmb code.                    *
       *                                                                *
       * The program is enrolled in MPI (for a parallel executable) and *
       * the command line arguments are extracted after which the       *
       * Fortran main is called to do the actual job. The command line  *
       * arguments is the reason why C is used for the main; in Fortran *
       * command line arguments are an extension to the standard and    *
       * every MPI handles them differently. In C there is a unique     *
       * standard.                                                      *
       *                                                                *
       * Usage sequential mode:                                         *
       *       SUmb <parameter file>                                    *
       *                                                                *
       * Usage parallel mode:                                           *
       *       mpirun -np <nproc> SUmb <parameter file>                 *
       *                                                                *
       * Here parameter file is the name of the file with the input     *
       * parameters.                                                    *
       *                                                                *
       ******************************************************************
*/

int main(int argc, char **argv)
{
  int  len_exec, len_param;
  int  size_int;
  char *param_file;

  /* Enroll the program in MPI; for a sequential code this is a dummy. */

  MPI_Init(&argc, &argv);

  /* Put the command line arguments in a usuable form. Note that in C */
  /* the first argument is the name of the executable. Command line   */
  /* arguments beyond the color are ignored.                          */

  /* Name of the executable. */

  if(argc > 3) argc = 3;
  len_exec = strlen(argv[0]);

  /* The parameter file. If nothing was specified have it point to a */
  /* dummy string to avoid any problems when calling the Fortran.    */

  if(argc > 1)
  {
    param_file = argv[1];
    len_param  = strlen(param_file);
  }
  else
  {
    param_file = " ";
    len_param  = 1;
  }

  /* Determine the size of the integer and real types */
  /* used in the C part of this executable.           */
  /* This should be compatible with the Fortran part. */

  size_int = sizeof(SUmb_intT);

  /* Call the main Fortran routine to do the work. */

#ifdef FORTRAN_CAPITALS

  SUMB_FORTRAN(argv[0], &len_exec, param_file, &len_param,
	       &argc, &size_int, len_exec, len_param);

#elif FORTRAN_DOUBLE_UNDERSCORE

  sumb_fortran__(argv[0], &len_exec, param_file, &len_param,
		 &argc, &size_int, len_exec, len_param);

#elif FORTRAN_NO_UNDERSCORE
  
  sumb_fortran(argv[0], &len_exec, param_file, &len_param,
               &argc, &size_int, len_exec, len_param);

#else

  sumb_fortran_(argv[0], &len_exec, param_file, &len_param,
		&argc, &size_int, len_exec, len_param);

#endif

  /* Quit from MPI; for a sequential code this is a dummy. */

  MPI_Finalize();

  /* Return 0 to indicate succes. */

  return 0;
}
